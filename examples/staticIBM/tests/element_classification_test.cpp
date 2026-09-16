#include <cstdint>

using dfloat = double;
using dlong = std::int64_t;

#include "../../../src/ibm/elementClassification.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{
ibm::Triangle triangle(const std::array<double, 3> & a,
                       const std::array<double, 3> & b,
                       const std::array<double, 3> & c,
                       const std::size_t source)
{
  ibm::Triangle result;
  result.vertex[0] = {a[0], a[1], a[2]};
  result.vertex[1] = {b[0], b[1], b[2]};
  result.vertex[2] = {c[0], c[1], c[2]};
  result.sourceTriangle = source;
  return result;
}

std::vector<ibm::Triangle> cubeTriangles()
{
  const std::array<double, 3> p000 = {0.0, 0.0, 0.0};
  const std::array<double, 3> p100 = {1.0, 0.0, 0.0};
  const std::array<double, 3> p010 = {0.0, 1.0, 0.0};
  const std::array<double, 3> p110 = {1.0, 1.0, 0.0};
  const std::array<double, 3> p001 = {0.0, 0.0, 1.0};
  const std::array<double, 3> p101 = {1.0, 0.0, 1.0};
  const std::array<double, 3> p011 = {0.0, 1.0, 1.0};
  const std::array<double, 3> p111 = {1.0, 1.0, 1.0};

  std::vector<ibm::Triangle> triangles;
  auto add = [&](const std::array<double, 3> & a,
                 const std::array<double, 3> & b,
                 const std::array<double, 3> & c) {
    triangles.push_back(triangle(a, b, c, triangles.size()));
  };

  add(p000, p010, p110);
  add(p000, p110, p100);
  add(p001, p101, p111);
  add(p001, p111, p011);
  add(p000, p100, p101);
  add(p000, p101, p001);
  add(p010, p011, p111);
  add(p010, p111, p110);
  add(p000, p001, p011);
  add(p000, p011, p010);
  add(p100, p110, p111);
  add(p100, p111, p101);
  return triangles;
}

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

int main()
{
  const auto translation = ibm::parseTranslation("(1.0, -2.0, 0.5)");
  assert(translation[0] == 1.0);
  assert(translation[1] == -2.0);
  assert(translation[2] == 0.5);

  bool rejected = false;
  try
  {
    (void)ibm::parseTranslation("1 2");
  }
  catch (const std::invalid_argument &)
  {
    rejected = true;
  }
  assert(rejected);

  auto triangles = cubeTriangles();
  auto shifted = triangles;
  ibm::translateTriangles(shifted, translation);
  assert(std::abs(shifted.front().vertex[0].x - (triangles.front().vertex[0].x + 1.0)) <
         1e-15);
  assert(std::abs(shifted.front().vertex[0].y - (triangles.front().vertex[0].y - 2.0)) <
         1e-15);
  assert(std::abs(shifted.front().vertex[0].z - (triangles.front().vertex[0].z + 0.5)) <
         1e-15);

  std::vector<dfloat> x;
  std::vector<dfloat> y;
  std::vector<dfloat> z;
  appendBox({2.0, 2.0, 2.0}, {3.0, 3.0, 3.0}, x, y, z);
  appendBox({0.25, 0.25, 0.25}, {0.75, 0.75, 0.75}, x, y, z);
  appendBox({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, x, y, z);
  appendBox({1.0, 0.25, 0.25}, {1.5, 0.75, 0.75}, x, y, z);

  const auto regions = ibm::classifyElements(triangles, x, y, z, 4, 8);
  assert(regions[0] == ibm::ElementRegion::Fluid);
  assert(regions[1] == ibm::ElementRegion::Solid);
  assert(regions[2] == ibm::ElementRegion::Cut);
  assert(regions[3] == ibm::ElementRegion::Cut);

  const auto counts = ibm::countElementRegions(regions);
  assert(counts.fluid == 1);
  assert(counts.solid == 1);
  assert(counts.cut == 2);

  const auto field = ibm::makeElementRegionField(regions, 8);
  assert(field.size() == 32);
  assert(field[0] == 0.0);
  assert(field[8] == 2.0);
  assert(field[16] == 1.0);
  assert(field[24] == 1.0);
  return 0;
}
