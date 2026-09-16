#include <cstdint>

using dfloat = double;
using dlong = std::int64_t;

#include "../../../src/ibm/elementClassification.hpp"
#include "../../../src/ibm/stl_reader.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace
{
void appendBox(const std::array<double, 3> & minimum,
               const std::array<double, 3> & maximum,
               std::vector<dfloat> & x,
               std::vector<dfloat> & y,
               std::vector<dfloat> & z)
{
  const std::array<std::array<double, 3>, 8> nodes = {{
      {minimum[0], minimum[1], minimum[2]},
      {maximum[0], minimum[1], minimum[2]},
      {maximum[0], maximum[1], minimum[2]},
      {minimum[0], maximum[1], minimum[2]},
      {minimum[0], minimum[1], maximum[2]},
      {maximum[0], minimum[1], maximum[2]},
      {maximum[0], maximum[1], maximum[2]},
      {minimum[0], maximum[1], maximum[2]},
  }};
  for (const auto & node : nodes)
  {
    x.push_back(node[0]);
    y.push_back(node[1]);
    z.push_back(node[2]);
  }
}
}

int main(int argc, char ** argv)
{
  assert(argc == 2);
  auto triangles = ibm::readBinaryStl(argv[1]);
  assert(triangles.size() == 166684);

  const auto translation = ibm::parseTranslation("0.675, 0.675, 0.7575");
  ibm::translateTriangles(triangles, translation);

  std::array<double, 3> minimum = {
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max()};
  std::array<double, 3> maximum = {
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::lowest()};
  for (const auto & triangle : triangles)
    for (const auto & vertex : triangle.vertex)
    {
      minimum[0] = std::min(minimum[0], vertex.x);
      minimum[1] = std::min(minimum[1], vertex.y);
      minimum[2] = std::min(minimum[2], vertex.z);
      maximum[0] = std::max(maximum[0], vertex.x);
      maximum[1] = std::max(maximum[1], vertex.y);
      maximum[2] = std::max(maximum[2], vertex.z);
    }

  assert(minimum[0] > 0.0);
  assert(minimum[1] > 0.0);
  assert(minimum[2] > 0.0);
  assert(maximum[0] < 1.35);
  assert(maximum[1] < 1.35);
  assert(maximum[2] < 2.70);

  // A finer box than the runnable example is used here so the gyroid has
  // elements that are entirely inside the closed STL volume.
  constexpr int nx = 32;
  constexpr int ny = 32;
  constexpr int nz = 64;
  std::vector<dfloat> x;
  std::vector<dfloat> y;
  std::vector<dfloat> z;
  x.reserve(nx * ny * nz * 8);
  y.reserve(nx * ny * nz * 8);
  z.reserve(nx * ny * nz * 8);
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
      {
        const double x0 = 1.35 * i / nx;
        const double x1 = 1.35 * (i + 1) / nx;
        const double y0 = 1.35 * j / ny;
        const double y1 = 1.35 * (j + 1) / ny;
        const double z0 = 2.70 * k / nz;
        const double z1 = 2.70 * (k + 1) / nz;
        appendBox({x0, y0, z0}, {x1, y1, z1}, x, y, z);
      }

  const auto regions = ibm::classifyElements(triangles, x, y, z, nx * ny * nz, 8);
  const auto counts = ibm::countElementRegions(regions);
  assert(counts.fluid > 0);
  assert(counts.cut > 0);
  assert(counts.solid > 0);
  return 0;
}
