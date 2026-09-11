#include <cstdint>

#ifdef TEST_DFLOAT_FLOAT
using dfloat = float;
#else
using dfloat = double;
#endif
using dlong = std::int64_t;

#include "../ibmGeometry.hpp"
#include "../stl_reader.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <type_traits>

int main(int argc, char ** argv)
{
  assert(argc == 2);
  const auto triangles = ibm::readBinaryStl(argv[1]);
  assert(triangles.size() == 4);

  ibm::SamplingStatistics coarseStatistics;
  const auto coarse = ibm::sampleStaticSurface(triangles, 0.40, 0.05, 30,
                                                &coarseStatistics);
  ibm::validateMarkers(triangles, coarse, 0.05);
  assert(coarse.size() == 16);
  assert(coarseStatistics.maximumLeafMaxEdge <= 0.40);

  ibm::SamplingStatistics refinedStatistics;
  const auto refined = ibm::sampleStaticSurface(triangles, 0.20, 0.05, 30,
                                                 &refinedStatistics);
  ibm::validateMarkers(triangles, refined, 0.05);
  assert(refined.size() == 64);
  assert(refined.size() > coarse.size());
  assert(refinedStatistics.maximumLeafMaxEdge <= 0.20);

  double area = 0.0;
  for (const auto & triangle : triangles)
    area += ibm::triangleArea(triangle);
  double volume = 0.0;
  for (const auto value : refined.volume)
    volume += value;
  const double tolerance = std::is_same<dfloat, float>::value ? 1e-6 : 1e-13;
  assert(std::abs(volume - area * 0.05) < tolerance);

  bool rejected = false;
  auto degenerate = triangles;
  degenerate[0].vertex[2] = degenerate[0].vertex[1];
  try
  {
    (void)ibm::sampleStaticSurface(degenerate, 0.40, 0.05);
  }
  catch (const std::runtime_error &)
  {
    rejected = true;
  }
  assert(rejected);
  return 0;
}
