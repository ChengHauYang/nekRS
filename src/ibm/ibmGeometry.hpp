#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace ibm
{

struct Vec3
{
  double x;
  double y;
  double z;
};

struct Triangle
{
  Vec3 vertex[3];
  std::size_t sourceTriangle;
};

struct MarkerSet
{
  // Component-major: [x0 ... xN-1 | y0 ... yN-1 | z0 ... zN-1].
  std::vector<dfloat> coordinates;

  // Component-major no-slip target velocity.
  std::vector<dfloat> targetVelocity;

  // Marker volume: dVp = leaf triangle area * shellThickness.
  std::vector<dfloat> volume;

  // Original STL triangle that produced each marker.
  std::vector<dlong> sourceTriangle;

  std::size_t size() const { return volume.size(); }
};

struct SamplingStatistics
{
  double minimumLeafMaxEdge = std::numeric_limits<double>::infinity();
  double maximumLeafMaxEdge = 0.0;
};

namespace detail
{

inline Vec3 subtract(const Vec3 & a, const Vec3 & b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

inline Vec3 midpoint(const Vec3 & a, const Vec3 & b)
{
  return {0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z)};
}

inline Vec3 cross(const Vec3 & a, const Vec3 & b)
{
  return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

inline double norm(const Vec3 & value)
{
  return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

inline double edgeLength(const Vec3 & a, const Vec3 & b)
{
  return norm(subtract(a, b));
}

inline double triangleArea(const Triangle & triangle)
{
  return 0.5 * norm(cross(subtract(triangle.vertex[1], triangle.vertex[0]),
                          subtract(triangle.vertex[2], triangle.vertex[0])));
}

inline double maxEdgeLength(const Triangle & triangle)
{
  return std::max({edgeLength(triangle.vertex[0], triangle.vertex[1]),
                   edgeLength(triangle.vertex[1], triangle.vertex[2]),
                   edgeLength(triangle.vertex[2], triangle.vertex[0])});
}

inline Vec3 centroid(const Triangle & triangle)
{
  return {(triangle.vertex[0].x + triangle.vertex[1].x + triangle.vertex[2].x) / 3.0,
          (triangle.vertex[0].y + triangle.vertex[1].y + triangle.vertex[2].y) / 3.0,
          (triangle.vertex[0].z + triangle.vertex[1].z + triangle.vertex[2].z) / 3.0};
}

inline bool finite(const Vec3 & value)
{
  return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

inline void sampleTriangle(const Triangle & triangle,
                           const double originalMaxEdge,
                           const double markerSpacing,
                           const double shellThickness,
                           const unsigned int depth,
                           const unsigned int maxSubdivisionDepth,
                           std::vector<Vec3> & positions,
                           std::vector<double> & volumes,
                           std::vector<std::size_t> & sources,
                           SamplingStatistics * statistics)
{
  const double currentMaxEdge = maxEdgeLength(triangle);
  if (currentMaxEdge <= markerSpacing)
  {
    const double area = triangleArea(triangle);
    const double volume = area * shellThickness;
    const Vec3 position = centroid(triangle);
    if (!finite(position) || !std::isfinite(volume) || volume <= 0.0)
      throw std::runtime_error("IBM sampling generated a non-finite or non-positive marker");

    positions.push_back(position);
    volumes.push_back(volume);
    sources.push_back(triangle.sourceTriangle);
    if (statistics)
    {
      statistics->minimumLeafMaxEdge =
          std::min(statistics->minimumLeafMaxEdge, currentMaxEdge);
      statistics->maximumLeafMaxEdge =
          std::max(statistics->maximumLeafMaxEdge, currentMaxEdge);
    }
    return;
  }

  if (depth >= maxSubdivisionDepth)
  {
    std::ostringstream message;
    message << "IBM subdivision depth limit reached for source triangle "
            << triangle.sourceTriangle << ": original max edge=" << originalMaxEdge
            << ", current max edge=" << currentMaxEdge
            << ", marker spacing=" << markerSpacing << ", depth=" << depth;
    throw std::runtime_error(message.str());
  }

  const Vec3 & v0 = triangle.vertex[0];
  const Vec3 & v1 = triangle.vertex[1];
  const Vec3 & v2 = triangle.vertex[2];
  const Vec3 m01 = midpoint(v0, v1);
  const Vec3 m12 = midpoint(v1, v2);
  const Vec3 m20 = midpoint(v2, v0);

  // Fixed child order and recursive DFS make marker ordering deterministic.
  const std::array<Triangle, 4> children = {{
      {{v0, m01, m20}, triangle.sourceTriangle},
      {{m01, v1, m12}, triangle.sourceTriangle},
      {{m20, m12, v2}, triangle.sourceTriangle},
      {{m01, m12, m20}, triangle.sourceTriangle},
  }};

  for (const auto & child : children)
    sampleTriangle(child,
                   originalMaxEdge,
                   markerSpacing,
                   shellThickness,
                   depth + 1,
                   maxSubdivisionDepth,
                   positions,
                   volumes,
                   sources,
                   statistics);
}

} // namespace detail

inline double triangleArea(const Triangle & triangle)
{
  return detail::triangleArea(triangle);
}

inline MarkerSet sampleStaticSurface(const std::vector<Triangle> & triangles,
                                     const double markerSpacing,
                                     const double shellThickness,
                                     const unsigned int maxSubdivisionDepth = 30,
                                     SamplingStatistics * statistics = nullptr)
{
  if (!std::isfinite(markerSpacing) || markerSpacing <= 0.0)
    throw std::invalid_argument("markerSpacing must be finite and positive");
  if (!std::isfinite(shellThickness) || shellThickness <= 0.0)
    throw std::invalid_argument("shellThickness must be finite and positive");
  if (triangles.empty())
    throw std::invalid_argument("cannot sample an empty STL surface");

  std::vector<Vec3> positions;
  std::vector<double> volumes;
  std::vector<std::size_t> sources;

  for (const auto & triangle : triangles)
  {
    if (triangle.sourceTriangle >= triangles.size())
      throw std::runtime_error("source triangle index is out of range");
    const double area = detail::triangleArea(triangle);
    const double originalMaxEdge = detail::maxEdgeLength(triangle);
    if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(originalMaxEdge))
      throw std::runtime_error("input contains an invalid triangle");
    detail::sampleTriangle(triangle,
                           originalMaxEdge,
                           markerSpacing,
                           shellThickness,
                           0,
                           maxSubdivisionDepth,
                           positions,
                           volumes,
                           sources,
                           statistics);
  }

  const std::size_t count = positions.size();
  if (count > static_cast<std::size_t>(std::numeric_limits<dlong>::max()))
    throw std::overflow_error("marker count does not fit in dlong");
  if (count > std::numeric_limits<std::size_t>::max() / 3)
    throw std::overflow_error("component-major marker array size overflow");

  MarkerSet markers;
  markers.coordinates.resize(3 * count);
  markers.targetVelocity.assign(3 * count, static_cast<dfloat>(0));
  markers.volume.resize(count);
  markers.sourceTriangle.resize(count);

  for (std::size_t marker = 0; marker < count; ++marker)
  {
    markers.coordinates[marker] = static_cast<dfloat>(positions[marker].x);
    markers.coordinates[marker + count] = static_cast<dfloat>(positions[marker].y);
    markers.coordinates[marker + 2 * count] = static_cast<dfloat>(positions[marker].z);
    markers.volume[marker] = static_cast<dfloat>(volumes[marker]);
    markers.sourceTriangle[marker] = static_cast<dlong>(sources[marker]);

    // Catch values that were finite in double but overflowed or underflowed in dfloat.
    if (!std::isfinite(markers.coordinates[marker]) ||
        !std::isfinite(markers.coordinates[marker + count]) ||
        !std::isfinite(markers.coordinates[marker + 2 * count]) ||
        !std::isfinite(markers.volume[marker]) || markers.volume[marker] <= 0)
      throw std::runtime_error("marker cannot be represented safely with dfloat");
  }
  return markers;
}

inline void validateMarkers(const std::vector<Triangle> & triangles,
                            const MarkerSet & markers,
                            const double shellThickness)
{
  const std::size_t count = markers.size();
  if (!std::isfinite(shellThickness) || shellThickness <= 0.0)
    throw std::invalid_argument("shellThickness must be finite and positive");
  if (count == 0 || markers.coordinates.size() != 3 * count ||
      markers.targetVelocity.size() != 3 * count ||
      markers.sourceTriangle.size() != count)
    throw std::runtime_error("MarkerSet arrays have inconsistent sizes");

  std::vector<double> sampledArea(triangles.size(), 0.0);
  for (std::size_t marker = 0; marker < count; ++marker)
  {
    const dlong source = markers.sourceTriangle[marker];
    if (source < 0 || static_cast<std::size_t>(source) >= triangles.size())
      throw std::runtime_error("marker source triangle index is out of range");
    if (!std::isfinite(markers.volume[marker]) || markers.volume[marker] <= 0)
      throw std::runtime_error("marker volume must be finite and positive");
    for (unsigned int component = 0; component < 3; ++component)
    {
      if (!std::isfinite(markers.coordinates[marker + component * count]) ||
          !std::isfinite(markers.targetVelocity[marker + component * count]))
        throw std::runtime_error("marker coordinates and targets must be finite");
      if (markers.targetVelocity[marker + component * count] != 0)
        throw std::runtime_error("static no-slip target velocity must be zero");
    }
    sampledArea[static_cast<std::size_t>(source)] +=
        static_cast<double>(markers.volume[marker]) / shellThickness;
  }

  double originalTotalArea = 0.0;
  double sampledTotalArea = 0.0;
  for (std::size_t source = 0; source < triangles.size(); ++source)
  {
    const double originalArea = detail::triangleArea(triangles[source]);
    const double tolerance = 256.0 * std::numeric_limits<dfloat>::epsilon() *
                             std::max(1.0, originalArea);
    if (std::abs(sampledArea[source] - originalArea) > tolerance)
    {
      std::ostringstream message;
      message << "marker area conservation failed for source triangle " << source
              << ": sampled=" << sampledArea[source] << ", original=" << originalArea
              << ", tolerance=" << tolerance;
      throw std::runtime_error(message.str());
    }
    originalTotalArea += originalArea;
    sampledTotalArea += sampledArea[source];
  }

  const double globalTolerance = 256.0 * std::numeric_limits<dfloat>::epsilon() *
                                 std::max(1.0, originalTotalArea);
  if (std::abs(sampledTotalArea - originalTotalArea) > globalTolerance)
    throw std::runtime_error("global marker area/volume conservation failed");
}

} // namespace ibm
