#pragma once

#include "ibmGeometry.hpp"
#include <array>
#include <cctype>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>

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

struct Bin3D
{
  long long x;
  long long y;
  long long z;

  bool operator==(const Bin3D & other) const
  {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct BinHasher
{
  std::size_t operator()(const Bin3D & bin) const
  {
    const std::size_t x = std::hash<long long>{}(bin.x);
    const std::size_t y = std::hash<long long>{}(bin.y);
    const std::size_t z = std::hash<long long>{}(bin.z);
    return x ^ (y << 1) ^ (z << 2);
  }
};

struct PeriodicDomain
{
  std::array<bool, 3> enabled = {false, false, false};
  std::array<double, 3> minimum = {0.0, 0.0, 0.0};
  std::array<double, 3> length = {0.0, 0.0, 0.0};
};

struct ImageNode
{
  std::size_t index;
  double x;
  double y;
  double z;
};

inline PeriodicDomain periodicDomain(const std::vector<dfloat> & node_x,
                                     const std::vector<dfloat> & node_y,
                                     const std::vector<dfloat> & node_z,
                                     const double supportRadius,
                                     const std::string & periodicDimensions,
                                     MPI_Comm comm)
{
  PeriodicDomain domain;
  std::string dimensions;
  dimensions.reserve(periodicDimensions.size());
  for (const unsigned char value : periodicDimensions)
    if (!std::isspace(value) && value != ',')
      dimensions.push_back(static_cast<char>(std::tolower(value)));

  if (dimensions.empty() || dimensions == "none")
    return domain;

  for (const char dimension : dimensions)
  {
    if (dimension < 'x' || dimension > 'z')
      throw std::invalid_argument(
          "periodic_dimensions must be none or a combination of x, y, and z");
    domain.enabled[dimension - 'x'] = true;
  }

  std::array<double, 3> localMinimum = {
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()};
  std::array<double, 3> localMaximum = {
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()};
  for (std::size_t n = 0; n < node_x.size(); ++n)
  {
    localMinimum[0] = std::min(localMinimum[0], static_cast<double>(node_x[n]));
    localMinimum[1] = std::min(localMinimum[1], static_cast<double>(node_y[n]));
    localMinimum[2] = std::min(localMinimum[2], static_cast<double>(node_z[n]));
    localMaximum[0] = std::max(localMaximum[0], static_cast<double>(node_x[n]));
    localMaximum[1] = std::max(localMaximum[1], static_cast<double>(node_y[n]));
    localMaximum[2] = std::max(localMaximum[2], static_cast<double>(node_z[n]));
  }

  std::array<double, 3> globalMaximum;
  MPI_Allreduce(localMinimum.data(), domain.minimum.data(), 3, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(localMaximum.data(), globalMaximum.data(), 3, MPI_DOUBLE, MPI_MAX, comm);

  for (int component = 0; component < 3; ++component)
  {
    if (!domain.enabled[component])
      continue;
    domain.length[component] = globalMaximum[component] - domain.minimum[component];
    if (!std::isfinite(domain.length[component]) ||
        domain.length[component] <= 2.0 * supportRadius)
    {
      std::ostringstream message;
      message << "periodic domain length in " << static_cast<char>('x' + component)
              << " must be greater than twice the IBM support radius: length="
              << domain.length[component] << ", R=" << supportRadius;
      throw std::invalid_argument(message.str());
    }
  }
  return domain;
}

} // namespace detail

inline InteractionMaps buildInteractionMaps(
    const MarkerSet& markers,
    const std::vector<dfloat>& node_x,
    const std::vector<dfloat>& node_y,
    const std::vector<dfloat>& node_z,
    const std::vector<dfloat>& node_jw,
    const double R,
    const double sigma,
    MPI_Comm comm,
    const std::string & periodicDimensions = "none")
{
  if (R <= 0.0 || sigma <= 0.0) {
    throw std::invalid_argument("R and sigma must be positive");
  }

  const double R2 = R * R;
  const double sigma2 = sigma * sigma;
  const std::size_t Np = markers.size();
  const std::size_t Nlocal = node_x.size();
  if (node_y.size() != Nlocal || node_z.size() != Nlocal || node_jw.size() != Nlocal)
    throw std::invalid_argument("IBM node coordinate and quadrature arrays must have equal sizes");

  const auto periodic =
      detail::periodicDomain(node_x, node_y, node_z, R, periodicDimensions, comm);

  std::unordered_map<detail::Bin3D,
                     std::vector<detail::ImageNode>,
                     detail::BinHasher> bins;
  for (std::size_t n = 0; n < Nlocal; ++n)
  {
    for (int imageX = periodic.enabled[0] ? -1 : 0;
         imageX <= (periodic.enabled[0] ? 1 : 0);
         ++imageX)
      for (int imageY = periodic.enabled[1] ? -1 : 0;
           imageY <= (periodic.enabled[1] ? 1 : 0);
           ++imageY)
        for (int imageZ = periodic.enabled[2] ? -1 : 0;
             imageZ <= (periodic.enabled[2] ? 1 : 0);
             ++imageZ)
        {
          const detail::ImageNode image = {
              n,
              node_x[n] + imageX * periodic.length[0],
              node_y[n] + imageY * periodic.length[1],
              node_z[n] + imageZ * periodic.length[2]};
          const detail::Bin3D bin = {
              static_cast<long long>(std::floor(image.x / R)),
              static_cast<long long>(std::floor(image.y / R)),
              static_cast<long long>(std::floor(image.z / R))};
          bins[bin].push_back(image);
        }
  }

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

    const detail::Bin3D markerBin = {
        static_cast<long long>(std::floor(px / R)),
        static_cast<long long>(std::floor(py / R)),
        static_cast<long long>(std::floor(pz / R))};

    double local_sp = 0.0;
    auto& conn = particle_connections[p];

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const detail::Bin3D searchBin = {
              markerBin.x + dx, markerBin.y + dy, markerBin.z + dz};
          const auto bin = bins.find(searchBin);
          if (bin == bins.end())
            continue;

          for (const auto & image : bin->second) {
            const std::size_t n = image.index;
            const double dist2 = (px - image.x) * (px - image.x) +
                                 (py - image.y) * (py - image.y) +
                                 (pz - image.z) * (pz - image.z);
            if (dist2 <= R2) {
              const double G = detail::gaussian(dist2, sigma2);
              local_sp += G * node_jw[n];
              conn.push_back({n, G, node_jw[n]});
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
