#pragma once

#include "ibmGeometry.hpp"
#include <vector>
#include <cmath>
#include <mpi.h>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <numeric>

namespace ibm
{

struct CsrMap
{
  std::vector<dlong> offsets;
  std::vector<dlong> indices;
  std::vector<dfloat> weights;
};

struct InteractionMaps
{
  CsrMap interp;
  CsrMap spread;
};

namespace detail
{

inline double gaussian(double r2, double sigma2)
{
  return std::exp(-r2 / (2.0 * sigma2)) / std::pow(2.0 * M_PI * sigma2, 1.5);
}

struct Bin3D {
  long long x, y, z;
  bool operator==(const Bin3D& other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct BinHasher {
  std::size_t operator()(const Bin3D& b) const {
    std::size_t h1 = std::hash<long long>{}(b.x);
    std::size_t h2 = std::hash<long long>{}(b.y);
    std::size_t h3 = std::hash<long long>{}(b.z);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

} // namespace detail

inline InteractionMaps buildInteractionMaps(
    const MarkerSet& markers,
    const std::vector<dfloat>& node_x,
    const std::vector<dfloat>& node_y,
    const std::vector<dfloat>& node_z,
    const std::vector<dfloat>& node_jw,
    const double R,
    const double sigma,
    MPI_Comm comm)
{
  if (R <= 0.0 || sigma <= 0.0) {
    throw std::invalid_argument("R and sigma must be positive");
  }

  const double R2 = R * R;
  const double sigma2 = sigma * sigma;
  const std::size_t Np = markers.size();
  const std::size_t Nlocal = node_x.size();

  std::unordered_map<detail::Bin3D, std::vector<std::size_t>, detail::BinHasher> bins;
  for (std::size_t n = 0; n < Nlocal; ++n) {
    detail::Bin3D bin = {
      static_cast<long long>(std::floor(node_x[n] / R)),
      static_cast<long long>(std::floor(node_y[n] / R)),
      static_cast<long long>(std::floor(node_z[n] / R))
    };
    bins[bin].push_back(n);
  }

  std::vector<std::size_t> row_counts(Np, 0);
  std::vector<double> s_p_local(Np, 0.0);

  struct Connection {
    std::size_t n;
    double G;
    double dV_n;
  };
  std::vector<std::vector<Connection>> particle_connections(Np);

  for (std::size_t p = 0; p < Np; ++p) {
    double px = markers.coordinates[p];
    double py = markers.coordinates[p + Np];
    double pz = markers.coordinates[p + 2 * Np];

    detail::Bin3D pbin = {
      static_cast<long long>(std::floor(px / R)),
      static_cast<long long>(std::floor(py / R)),
      static_cast<long long>(std::floor(pz / R))
    };

    double local_sp = 0.0;
    auto& conn = particle_connections[p];

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          detail::Bin3D search_bin = {pbin.x + dx, pbin.y + dy, pbin.z + dz};
          auto it = bins.find(search_bin);
          if (it != bins.end()) {
            for (std::size_t n : it->second) {
              double dist2 = (px - node_x[n])*(px - node_x[n]) + 
                             (py - node_y[n])*(py - node_y[n]) + 
                             (pz - node_z[n])*(pz - node_z[n]);
              if (dist2 <= R2) {
                double G = detail::gaussian(dist2, sigma2);
                local_sp += G * node_jw[n];
                conn.push_back({n, G, node_jw[n]});
              }
            }
          }
        }
      }
    }
    s_p_local[p] = local_sp;
  }

  std::vector<double> s_p_global(Np, 0.0);
  MPI_Allreduce(s_p_local.data(), s_p_global.data(), Np, MPI_DOUBLE, MPI_SUM, comm);

  const double support_tolerance = 1e-12;
  for (std::size_t p = 0; p < Np; ++p) {
    if (s_p_global[p] < support_tolerance) {
      throw std::runtime_error("Marker unsupported or below tolerance");
    }
  }

  InteractionMaps maps;
  maps.interp.offsets.assign(Np + 1, 0);
  std::size_t total_interp_entries = 0;
  for (std::size_t p = 0; p < Np; ++p) {
    total_interp_entries += particle_connections[p].size();
    maps.interp.offsets[p + 1] = total_interp_entries;
  }
  
  maps.interp.indices.reserve(total_interp_entries);
  maps.interp.weights.reserve(total_interp_entries);

  std::vector<std::size_t> spread_counts(Nlocal, 0);

  for (std::size_t p = 0; p < Np; ++p) {
    double sp = s_p_global[p];
    for (const auto& conn : particle_connections[p]) {
      maps.interp.indices.push_back(static_cast<dlong>(conn.n));
      maps.interp.weights.push_back(static_cast<dfloat>(conn.G * conn.dV_n / sp));
      spread_counts[conn.n]++;
    }
  }

  maps.spread.offsets.assign(Nlocal + 1, 0);
  for (std::size_t n = 0; n < Nlocal; ++n) {
    maps.spread.offsets[n + 1] = maps.spread.offsets[n] + spread_counts[n];
  }

  std::size_t total_spread_entries = maps.spread.offsets.back();
  maps.spread.indices.resize(total_spread_entries);
  maps.spread.weights.resize(total_spread_entries);

  std::vector<std::size_t> current_spread_offset(maps.spread.offsets.begin(), maps.spread.offsets.end() - 1);
  for (std::size_t p = 0; p < Np; ++p) {
    double sp = s_p_global[p];
    double dVp = markers.volume[p];
    for (const auto& conn : particle_connections[p]) {
      std::size_t pos = current_spread_offset[conn.n]++;
      maps.spread.indices[pos] = static_cast<dlong>(p);
      maps.spread.weights[pos] = static_cast<dfloat>(conn.G * dVp / sp);
    }
  }

  return maps;
}

} // namespace ibm
