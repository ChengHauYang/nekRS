#pragma once

#include "ibmGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace ibm
{

enum class ElementRegion : int
{
  Fluid = 0,
  Cut = 1,
  Solid = 2
};

struct ElementRegionCounts
{
  unsigned long long fluid = 0;
  unsigned long long cut = 0;
  unsigned long long solid = 0;
};

namespace classification_detail
{

struct Aabb
{
  Vec3 minimum = {std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity()};
  Vec3 maximum = {-std::numeric_limits<double>::infinity(),
                  -std::numeric_limits<double>::infinity(),
                  -std::numeric_limits<double>::infinity()};
};

inline double component(const Vec3 & value, const int axis)
{
  if (axis == 0)
    return value.x;
  if (axis == 1)
    return value.y;
  return value.z;
}

inline void expand(Aabb & box, const Vec3 & point)
{
  box.minimum.x = std::min(box.minimum.x, point.x);
  box.minimum.y = std::min(box.minimum.y, point.y);
  box.minimum.z = std::min(box.minimum.z, point.z);
  box.maximum.x = std::max(box.maximum.x, point.x);
  box.maximum.y = std::max(box.maximum.y, point.y);
  box.maximum.z = std::max(box.maximum.z, point.z);
}

inline void expand(Aabb & box, const Aabb & other)
{
  expand(box, other.minimum);
  expand(box, other.maximum);
}

inline Aabb triangleBounds(const Triangle & triangle)
{
  Aabb box;
  for (const auto & vertex : triangle.vertex)
    expand(box, vertex);
  return box;
}

inline Vec3 center(const Aabb & box)
{
  return {0.5 * (box.minimum.x + box.maximum.x),
          0.5 * (box.minimum.y + box.maximum.y),
          0.5 * (box.minimum.z + box.maximum.z)};
}

inline bool overlaps(const Aabb & first, const Aabb & second)
{
  const double scale =
      std::max({1.0,
                std::abs(first.minimum.x),
                std::abs(first.minimum.y),
                std::abs(first.minimum.z),
                std::abs(first.maximum.x),
                std::abs(first.maximum.y),
                std::abs(first.maximum.z),
                std::abs(second.minimum.x),
                std::abs(second.minimum.y),
                std::abs(second.minimum.z),
                std::abs(second.maximum.x),
                std::abs(second.maximum.y),
                std::abs(second.maximum.z)});
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  return first.maximum.x + tolerance >= second.minimum.x &&
         second.maximum.x + tolerance >= first.minimum.x &&
         first.maximum.y + tolerance >= second.minimum.y &&
         second.maximum.y + tolerance >= first.minimum.y &&
         first.maximum.z + tolerance >= second.minimum.z &&
         second.maximum.z + tolerance >= first.minimum.z;
}

inline double dot(const Vec3 & first, const Vec3 & second)
{
  return first.x * second.x + first.y * second.y + first.z * second.z;
}

inline bool separatedOnAxis(const std::array<Vec3, 3> & vertices,
                            const Vec3 & halfWidth,
                            const Vec3 & axis)
{
  const double axisNormSquared = dot(axis, axis);
  if (axisNormSquared <= std::numeric_limits<double>::epsilon())
    return false;

  const std::array<double, 3> projection = {
      dot(vertices[0], axis), dot(vertices[1], axis), dot(vertices[2], axis)};
  const auto range = std::minmax_element(projection.begin(), projection.end());
  const double radius = halfWidth.x * std::abs(axis.x) + halfWidth.y * std::abs(axis.y) +
                        halfWidth.z * std::abs(axis.z);
  const double scale =
      std::max({1.0, radius, std::abs(*range.first), std::abs(*range.second)});
  const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale;
  return *range.first > radius + tolerance || *range.second < -radius - tolerance;
}

// Separating-axis test between a triangle and an axis-aligned box.
inline bool triangleIntersectsBox(const Triangle & triangle, const Aabb & box)
{
  if (!overlaps(triangleBounds(triangle), box))
    return false;

  const Vec3 boxCenter = center(box);
  const Vec3 halfWidth = {0.5 * (box.maximum.x - box.minimum.x),
                          0.5 * (box.maximum.y - box.minimum.y),
                          0.5 * (box.maximum.z - box.minimum.z)};
  const std::array<Vec3, 3> vertices = {
      detail::subtract(triangle.vertex[0], boxCenter),
      detail::subtract(triangle.vertex[1], boxCenter),
      detail::subtract(triangle.vertex[2], boxCenter)};
  const std::array<Vec3, 3> edges = {detail::subtract(vertices[1], vertices[0]),
                                     detail::subtract(vertices[2], vertices[1]),
                                     detail::subtract(vertices[0], vertices[2])};
  const std::array<Vec3, 3> boxAxes = {{{1.0, 0.0, 0.0},
                                        {0.0, 1.0, 0.0},
                                        {0.0, 0.0, 1.0}}};

  for (const auto & axis : boxAxes)
    if (separatedOnAxis(vertices, halfWidth, axis))
      return false;

  const Vec3 triangleNormal = detail::cross(edges[0], edges[1]);
  if (separatedOnAxis(vertices, halfWidth, triangleNormal))
    return false;

  for (const auto & edge : edges)
    for (const auto & axis : boxAxes)
      if (separatedOnAxis(vertices, halfWidth, detail::cross(edge, axis)))
        return false;

  return true;
}

inline bool rayIntersectsBox(const Vec3 & origin, const Vec3 & direction, const Aabb & box)
{
  double lower = 0.0;
  double upper = std::numeric_limits<double>::infinity();
  for (int axis = 0; axis < 3; ++axis)
  {
    const double coordinate = component(origin, axis);
    const double slope = component(direction, axis);
    const double minimum = component(box.minimum, axis);
    const double maximum = component(box.maximum, axis);
    if (std::abs(slope) <= std::numeric_limits<double>::epsilon())
    {
      if (coordinate < minimum || coordinate > maximum)
        return false;
      continue;
    }

    double first = (minimum - coordinate) / slope;
    double second = (maximum - coordinate) / slope;
    if (first > second)
      std::swap(first, second);
    lower = std::max(lower, first);
    upper = std::min(upper, second);
    if (lower > upper)
      return false;
  }
  return upper > 0.0;
}

inline bool rayTriangleIntersection(const Vec3 & origin,
                                    const Vec3 & direction,
                                    const Triangle & triangle,
                                    double & distance)
{
  const Vec3 edge0 = detail::subtract(triangle.vertex[1], triangle.vertex[0]);
  const Vec3 edge1 = detail::subtract(triangle.vertex[2], triangle.vertex[0]);
  const Vec3 crossDirection = detail::cross(direction, edge1);
  const double determinant = dot(edge0, crossDirection);
  const double determinantTolerance =
      128.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, detail::norm(edge0) * detail::norm(edge1));
  if (std::abs(determinant) <= determinantTolerance)
    return false;

  const double inverseDeterminant = 1.0 / determinant;
  const Vec3 offset = detail::subtract(origin, triangle.vertex[0]);
  const double firstCoordinate = dot(offset, crossDirection) * inverseDeterminant;
  const double barycentricTolerance = 256.0 * std::numeric_limits<double>::epsilon();
  if (firstCoordinate < -barycentricTolerance || firstCoordinate > 1.0 + barycentricTolerance)
    return false;

  const Vec3 crossOffset = detail::cross(offset, edge0);
  const double secondCoordinate = dot(direction, crossOffset) * inverseDeterminant;
  if (secondCoordinate < -barycentricTolerance ||
      firstCoordinate + secondCoordinate > 1.0 + barycentricTolerance)
    return false;

  distance = dot(edge1, crossOffset) * inverseDeterminant;
  return distance > 256.0 * std::numeric_limits<double>::epsilon();
}

class TriangleTree
{
public:
  explicit TriangleTree(const std::vector<Triangle> & triangles) : _triangles(triangles)
  {
    if (_triangles.empty())
      throw std::invalid_argument("cannot classify elements with an empty triangle set");

    _triangleBounds.reserve(_triangles.size());
    _triangleCenters.reserve(_triangles.size());
    for (const auto & triangle : _triangles)
    {
      _triangleBounds.push_back(triangleBounds(triangle));
      _triangleCenters.push_back(center(_triangleBounds.back()));
    }
    _order.resize(_triangles.size());
    std::iota(_order.begin(), _order.end(), 0);
    _nodes.reserve(2 * _triangles.size());
    build(0, _order.size());
  }

  bool intersects(const Aabb & box) const { return intersects(0, box); }

  bool contains(const Vec3 & point) const
  {
    // Three non-axis-aligned rays avoid the common edge/vertex degeneracies of
    // parity ray casting. Coincident triangle hits are merged before counting.
    const std::array<Vec3, 3> directions = {{{1.0, 0.3713906763541037, 0.6947465906068658},
                                              {0.2386191860831969, 1.0, 0.5176380902050415},
                                              {0.6139406135149205, 0.2911437827766148, 1.0}}};
    int insideVotes = 0;
    for (const auto & direction : directions)
    {
      std::vector<double> hits;
      collectRayHits(0, point, direction, hits);
      std::sort(hits.begin(), hits.end());
      std::size_t uniqueHits = 0;
      double previous = 0.0;
      for (const double hit : hits)
      {
        const double tolerance =
            1024.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(hit));
        if (uniqueHits == 0 || std::abs(hit - previous) > tolerance)
        {
          ++uniqueHits;
          previous = hit;
        }
      }
      insideVotes += uniqueHits % 2;
    }
    return insideVotes >= 2;
  }

private:
  struct Node
  {
    Aabb bounds;
    std::size_t begin = 0;
    std::size_t end = 0;
    std::size_t left = 0;
    std::size_t right = 0;
    bool leaf = false;
  };

  std::size_t build(const std::size_t begin, const std::size_t end)
  {
    const std::size_t nodeIndex = _nodes.size();
    _nodes.emplace_back();
    Aabb bounds;
    Aabb centerBounds;
    for (std::size_t position = begin; position < end; ++position)
    {
      const std::size_t triangle = _order[position];
      expand(bounds, _triangleBounds[triangle]);
      expand(centerBounds, _triangleCenters[triangle]);
    }

    _nodes[nodeIndex].bounds = bounds;
    _nodes[nodeIndex].begin = begin;
    _nodes[nodeIndex].end = end;
    constexpr std::size_t leafSize = 8;
    if (end - begin <= leafSize)
    {
      _nodes[nodeIndex].leaf = true;
      return nodeIndex;
    }

    const std::array<double, 3> extent = {
        centerBounds.maximum.x - centerBounds.minimum.x,
        centerBounds.maximum.y - centerBounds.minimum.y,
        centerBounds.maximum.z - centerBounds.minimum.z};
    const int axis = static_cast<int>(
        std::distance(extent.begin(), std::max_element(extent.begin(), extent.end())));
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(_order.begin() + begin,
                     _order.begin() + middle,
                     _order.begin() + end,
                     [&](const std::size_t first, const std::size_t second) {
                       const double firstCoordinate = component(_triangleCenters[first], axis);
                       const double secondCoordinate = component(_triangleCenters[second], axis);
                       return firstCoordinate < secondCoordinate ||
                              (firstCoordinate == secondCoordinate && first < second);
                     });
    const std::size_t left = build(begin, middle);
    const std::size_t right = build(middle, end);
    _nodes[nodeIndex].left = left;
    _nodes[nodeIndex].right = right;
    return nodeIndex;
  }

  bool intersects(const std::size_t nodeIndex, const Aabb & box) const
  {
    const Node & node = _nodes[nodeIndex];
    if (!overlaps(node.bounds, box))
      return false;
    if (node.leaf)
    {
      for (std::size_t position = node.begin; position < node.end; ++position)
      {
        const std::size_t triangle = _order[position];
        if (triangleIntersectsBox(_triangles[triangle], box))
          return true;
      }
      return false;
    }
    return intersects(node.left, box) || intersects(node.right, box);
  }

  void collectRayHits(const std::size_t nodeIndex,
                      const Vec3 & origin,
                      const Vec3 & direction,
                      std::vector<double> & hits) const
  {
    const Node & node = _nodes[nodeIndex];
    if (!rayIntersectsBox(origin, direction, node.bounds))
      return;
    if (node.leaf)
    {
      for (std::size_t position = node.begin; position < node.end; ++position)
      {
        double distance = 0.0;
        if (rayTriangleIntersection(origin, direction, _triangles[_order[position]], distance))
          hits.push_back(distance);
      }
      return;
    }
    collectRayHits(node.left, origin, direction, hits);
    collectRayHits(node.right, origin, direction, hits);
  }

  const std::vector<Triangle> & _triangles;
  std::vector<Aabb> _triangleBounds;
  std::vector<Vec3> _triangleCenters;
  std::vector<std::size_t> _order;
  std::vector<Node> _nodes;
};

} // namespace classification_detail

inline std::vector<ElementRegion> classifyElements(const std::vector<Triangle> & triangles,
                                                   const std::vector<dfloat> & nodeX,
                                                   const std::vector<dfloat> & nodeY,
                                                   const std::vector<dfloat> & nodeZ,
                                                   const std::size_t elementCount,
                                                   const std::size_t nodesPerElement)
{
  const std::size_t nodeCount = elementCount * nodesPerElement;
  if (nodesPerElement == 0 || nodeX.size() != nodeCount || nodeY.size() != nodeCount ||
      nodeZ.size() != nodeCount)
    throw std::invalid_argument("invalid mesh arrays for IBM element classification");

  const classification_detail::TriangleTree tree(triangles);
  std::vector<ElementRegion> regions(elementCount, ElementRegion::Fluid);
  for (std::size_t element = 0; element < elementCount; ++element)
  {
    classification_detail::Aabb bounds;
    Vec3 centroid = {0.0, 0.0, 0.0};
    for (std::size_t localNode = 0; localNode < nodesPerElement; ++localNode)
    {
      const std::size_t node = element * nodesPerElement + localNode;
      const Vec3 point = {static_cast<double>(nodeX[node]),
                          static_cast<double>(nodeY[node]),
                          static_cast<double>(nodeZ[node])};
      classification_detail::expand(bounds, point);
      centroid.x += point.x;
      centroid.y += point.y;
      centroid.z += point.z;
    }
    centroid.x /= nodesPerElement;
    centroid.y /= nodesPerElement;
    centroid.z /= nodesPerElement;

    if (tree.intersects(bounds))
      regions[element] = ElementRegion::Cut;
    else if (tree.contains(centroid))
      regions[element] = ElementRegion::Solid;
  }
  return regions;
}

inline ElementRegionCounts countElementRegions(const std::vector<ElementRegion> & regions)
{
  ElementRegionCounts counts;
  for (const auto region : regions)
  {
    if (region == ElementRegion::Fluid)
      ++counts.fluid;
    else if (region == ElementRegion::Cut)
      ++counts.cut;
    else
      ++counts.solid;
  }
  return counts;
}

inline std::vector<dfloat> makeElementRegionField(const std::vector<ElementRegion> & regions,
                                                 const std::size_t nodesPerElement)
{
  std::vector<dfloat> field(regions.size() * nodesPerElement);
  for (std::size_t element = 0; element < regions.size(); ++element)
    std::fill(field.begin() + element * nodesPerElement,
              field.begin() + (element + 1) * nodesPerElement,
              static_cast<dfloat>(regions[element]));
  return field;
}

} // namespace ibm
