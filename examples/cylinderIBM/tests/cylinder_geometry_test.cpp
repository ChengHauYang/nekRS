#include <cstdint>

using dfloat = double;
using dlong = std::int64_t;

#include "../../../src/ibm/ibmGeometry.hpp"
#include "../../../src/ibm/stl_reader.hpp"

#include <cassert>
#include <cmath>

int main(int argc, char ** argv)
{
  assert(argc == 2);
  const auto triangles = ibm::readBinaryStl(argv[1]);
  assert(triangles.size() == 1024);

  double area = 0.0;
  for (const auto & triangle : triangles)
    area += ibm::triangleArea(triangle);
  const double exactArea = std::acos(-1.0);
  assert(std::abs(area - exactArea) / exactArea < 5e-4);

  constexpr double markerSpacing = 0.14;
  constexpr double shellThickness = 0.125;
  ibm::SamplingStatistics statistics;
  const auto markers = ibm::sampleStaticSurface(
      triangles, markerSpacing, shellThickness, 30, &statistics);
  ibm::validateMarkers(triangles, markers, shellThickness);
  assert(markers.size() == 1024);
  assert(statistics.maximumLeafMaxEdge <= markerSpacing);
  return 0;
}
