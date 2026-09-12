#include <cstdint>

using dfloat = double;
using dlong = std::int64_t;

#include "../../../src/ibm/ibmInteraction.hpp"

#include <algorithm>
#include <cassert>
#include <vector>

int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);

  ibm::MarkerSet markers;
  markers.coordinates = {0.0, 0.0, 0.98};
  markers.targetVelocity = {0.0, 0.0, 0.0};
  markers.volume = {1.0};
  markers.sourceTriangle = {0};

  const std::vector<dfloat> x(5, 0.0);
  const std::vector<dfloat> y(5, 0.0);
  const std::vector<dfloat> z = {0.0, 0.25, 0.50, 0.75, 1.0};
  const std::vector<dfloat> jw(5, 0.2);

  const auto nonperiodic =
      ibm::buildInteractionMaps(markers, x, y, z, jw, 0.1, 0.05, MPI_COMM_WORLD);
  assert(nonperiodic.interp.indices.size() == 1);
  assert(nonperiodic.interp.indices.front() == 4);

  const auto periodic =
      ibm::buildInteractionMaps(markers, x, y, z, jw, 0.1, 0.05, MPI_COMM_WORLD, "z");
  assert(periodic.interp.indices.size() == 2);
  assert(std::find(periodic.interp.indices.begin(), periodic.interp.indices.end(), 0) !=
         periodic.interp.indices.end());
  assert(std::find(periodic.interp.indices.begin(), periodic.interp.indices.end(), 4) !=
         periodic.interp.indices.end());

  MPI_Finalize();
  return 0;
}
