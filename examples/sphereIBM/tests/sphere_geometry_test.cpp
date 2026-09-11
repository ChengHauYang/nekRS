#include <cstdint>

using dfloat = double;
using dlong = std::int64_t;

#include "../../../src/ibm/ibmGeometry.hpp"
#include "../../../src/ibm/stl_reader.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

int main(int argc, char ** argv)
{
  assert(argc == 2);
  const auto triangles = ibm::readBinaryStl(argv[1]);
  assert(triangles.size() == 1520);

  std::array<double, 3> minimum = {
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max()};
  std::array<double, 3> maximum = {
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::lowest()};
  double area = 0;
  for (const auto & triangle : triangles)
  {
    area += ibm::triangleArea(triangle);
    for (const auto & vertex : triangle.vertex)
    {
      minimum[0] = std::min(minimum[0], vertex.x);
      minimum[1] = std::min(minimum[1], vertex.y);
      minimum[2] = std::min(minimum[2], vertex.z);
      maximum[0] = std::max(maximum[0], vertex.x);
      maximum[1] = std::max(maximum[1], vertex.y);
      maximum[2] = std::max(maximum[2], vertex.z);
    }
  }

  constexpr double coordinateTolerance = 1e-7;
  for (int component = 0; component < 3; ++component)
  {
    assert(std::abs(minimum[component] - 0.25) < coordinateTolerance);
    assert(std::abs(maximum[component] - 0.75) < coordinateTolerance);
  }

  const double exactArea = 4.0 * std::acos(-1.0) * 0.25 * 0.25;
  assert(std::abs(area - exactArea) / exactArea < 0.006);

  constexpr double markerSpacing = 0.04;
  constexpr double shellThickness = 0.04;
  ibm::SamplingStatistics statistics;
  const auto markers = ibm::sampleStaticSurface(
      triangles, markerSpacing, shellThickness, 30, &statistics);
  ibm::validateMarkers(triangles, markers, shellThickness);
  assert(markers.size() == 5840);
  assert(statistics.maximumLeafMaxEdge <= markerSpacing);
  return 0;
}
