#pragma once

#include "elementClassification.hpp"
#include "ibmGeometry.hpp"
#include "ibmInteraction.hpp"
#include "stl_reader.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ibm
{

struct StaticSurfaceDiagnostics
{
  dfloat maximumSlip = 0;
  dfloat rmsSlip = 0;
  dfloat force[3] = {0, 0, 0};
};

class StaticSurfaceIBM
{
public:
  void configure(setupAide & options)
  {
    platform->par->extract("casedata", "stl_file", _stlFile);
    platform->par->extract("casedata", "marker_spacing", _markerSpacing);
    platform->par->extract("casedata", "shell_thickness", _shellThickness);
    platform->par->extract("casedata", "support_radius", _supportRadius);
    platform->par->extract("casedata", "gaussian_width", _gaussianWidth);
    platform->par->extract("casedata", "max_subdivision_depth", _maxSubdivisionDepth);
    platform->par->extract("casedata", "periodic_dimensions", _periodicDimensions);
    std::string classifyElements;
    platform->par->extract("casedata", "classify_elements", classifyElements);
    if (!classifyElements.empty())
      _classifyElements = parseBool(classifyElements, "classify_elements");
    std::string stlTranslation;
    platform->par->extract("casedata", "stl_translation", stlTranslation);
    if (!stlTranslation.empty())
      _stlTranslation = parseTranslation(stlTranslation);
    options.getArgs("CI-MODE", _ciMode);

    if (_ciMode == 2)
      _markerSpacing *= 0.5;
  }

  void setup(nrs_t * nrs)
  {
    _nrs = nrs;
    const MPI_Comm comm = platform->comm.mpiComm();
    const int rank = platform->comm.mpiRank();
    std::string setupError;
    std::vector<Triangle> triangles;
    SamplingStatistics statistics;

    if (rank == 0)
    {
      try
      {
        validateParameters();
        triangles = readBinaryStl(_stlFile);
        translateTriangles(triangles, _stlTranslation);
        _markers = sampleStaticSurface(triangles,
                                       _markerSpacing,
                                       _shellThickness,
                                       static_cast<unsigned int>(_maxSubdivisionDepth),
                                       &statistics);
        validateMarkers(triangles, _markers, _shellThickness);
      }
      catch (const std::exception & error)
      {
        setupError = error.what();
      }
    }

    broadcastRootError(setupError, 0, comm);
    nekrsCheck(!setupError.empty(),
               comm,
               EXIT_FAILURE,
               "IBM geometry setup failed: %s\n",
               setupError.c_str());
    broadcastMarkerSet(_markers, 0, comm);
    validateBroadcastResult(_markers, comm);
    if (_classifyElements)
      broadcastTriangles(triangles, 0, comm);

    if (rank == 0)
      printGeometrySummary(triangles, statistics);

    buildDeviceData(triangles, comm);
  }

  void applyForcing(double)
  {
    if (_markers.size() == 0)
      return;

    interpolate(_nrs->fluid->o_U);

    // Direct forcing for the present first-order feedback prototype. Higher-order
    // BDF consistency is intentionally tracked separately from this refactor.
    const dfloat gainOverDt = 1.0 / _nrs->dt[0];
    ibm_compute_feedback(markerCount(),
                         gainOverDt,
                         _oMarkerVelocity,
                         _oTargetVelocity,
                         _oMarkerAcceleration);

    ibm_spread(_nrs->meshV->Nlocal,
               markerCount(),
               _nrs->fluid->fieldOffset,
               _oSpreadOffsets,
               _oSpreadIndices,
               _oSpreadWeights,
               _oMarkerAcceleration,
               _nrs->fluid->o_explicitTerms());
  }

  StaticSurfaceDiagnostics diagnostics()
  {
    interpolate(_nrs->fluid->o_U);
    _oMarkerVelocity.copyTo(_hostMarkerVelocity);

    StaticSurfaceDiagnostics result;
    double slipSquared = 0;
    for (std::size_t p = 0; p < _markers.size(); ++p)
    {
      double markerSlipSquared = 0;
      for (int component = 0; component < 3; ++component)
      {
        const std::size_t index = p + component * _markers.size();
        const double slip = _hostMarkerVelocity[index] - _markers.targetVelocity[index];
        markerSlipSquared += slip * slip;
      }
      result.maximumSlip = std::max(result.maximumSlip,
                                    static_cast<dfloat>(std::sqrt(markerSlipSquared)));
      slipSquared += markerSlipSquared;
    }
    result.rmsSlip = std::sqrt(slipSquared / _markers.size());

    std::vector<dfloat> acceleration(3 * _markers.size());
    _oMarkerAcceleration.copyTo(acceleration);
    for (std::size_t p = 0; p < _markers.size(); ++p)
      for (int component = 0; component < 3; ++component)
        result.force[component] +=
            acceleration[p + component * _markers.size()] * _markers.volume[p];
    return result;
  }

  const MarkerSet & markers() const { return _markers; }
  bool broadcastValidationPassed() const { return _broadcastValidationPassed; }
  bool hasElementClassification() const { return _classifyElements; }
  const deviceMemory<dfloat> & elementRegionField() const { return _oElementRegionField; }
  const ElementRegionCounts & elementRegionCounts() const { return _globalRegionCounts; }
  int ciMode() const { return _ciMode; }

private:
  dlong markerCount() const { return static_cast<dlong>(_markers.size()); }

  void validateParameters() const
  {
    if (_maxSubdivisionDepth < 0)
      throw std::invalid_argument("max_subdivision_depth must be non-negative");
    if (!std::isfinite(_supportRadius) || _supportRadius <= 0)
      throw std::invalid_argument("support_radius must be finite and positive");
    if (!std::isfinite(_gaussianWidth) || _gaussianWidth <= 0)
      throw std::invalid_argument("gaussian_width must be finite and positive");
  }

  static void checkMpiCount(const std::size_t count, const char * name, MPI_Comm comm)
  {
    nekrsCheck(count > static_cast<std::size_t>(INT_MAX),
               comm,
               EXIT_FAILURE,
               "IBM %s contains %zu entries, exceeding MPI int count\n",
               name,
               count);
  }

  static bool parseBool(const std::string & value, const char * name)
  {
    std::string cleaned;
    for (const char character : value)
      if (!std::isspace(static_cast<unsigned char>(character)))
        cleaned.push_back(static_cast<char>(
            std::tolower(static_cast<unsigned char>(character))));

    if (cleaned == "true" || cleaned == "1" || cleaned == "yes" ||
        cleaned == "on")
      return true;
    if (cleaned == "false" || cleaned == "0" || cleaned == "no" ||
        cleaned == "off")
      return false;

    std::ostringstream message;
    message << name << " must be true or false";
    throw std::invalid_argument(message.str());
  }

  static void broadcastRootError(std::string & message, const int root, MPI_Comm comm)
  {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    unsigned long long length =
        rank == root ? static_cast<unsigned long long>(message.size()) : 0;
    MPI_Bcast(&length, 1, MPI_UNSIGNED_LONG_LONG, root, comm);
    nekrsCheck(length > static_cast<unsigned long long>(INT_MAX),
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM setup error message exceeds MPI count\n");
    if (rank != root)
      message.resize(static_cast<std::size_t>(length));
    if (length)
      MPI_Bcast(&message[0], static_cast<int>(length), MPI_CHAR, root, comm);
  }

  static void broadcastMarkerSet(MarkerSet & values, const int root, MPI_Comm comm)
  {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    unsigned long long count =
        rank == root ? static_cast<unsigned long long>(values.size()) : 0;
    MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG_LONG, root, comm);

    nekrsCheck(count > static_cast<unsigned long long>(std::numeric_limits<dlong>::max()) ||
                   count > std::numeric_limits<std::size_t>::max() / 3,
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM marker count cannot be represented safely\n");
    const std::size_t markerCount = static_cast<std::size_t>(count);
    if (rank != root)
    {
      values.coordinates.resize(3 * markerCount);
      values.targetVelocity.resize(3 * markerCount);
      values.volume.resize(markerCount);
      values.sourceTriangle.resize(markerCount);
    }

    checkMpiCount(values.coordinates.size(), "coordinates", comm);
    checkMpiCount(values.targetVelocity.size(), "target velocities", comm);
    checkMpiCount(values.volume.size(), "volumes", comm);
    nekrsCheck(values.sourceTriangle.size() > static_cast<std::size_t>(INT_MAX) / sizeof(dlong),
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM source-triangle array exceeds MPI int byte count\n");

    MPI_Bcast(values.coordinates.data(),
              static_cast<int>(values.coordinates.size()),
              MPI_DFLOAT,
              root,
              comm);
    MPI_Bcast(values.targetVelocity.data(),
              static_cast<int>(values.targetVelocity.size()),
              MPI_DFLOAT,
              root,
              comm);
    MPI_Bcast(values.volume.data(),
              static_cast<int>(values.volume.size()),
              MPI_DFLOAT,
              root,
              comm);
    MPI_Bcast(values.sourceTriangle.data(),
              static_cast<int>(values.sourceTriangle.size() * sizeof(dlong)),
              MPI_BYTE,
              root,
              comm);
  }

  static void broadcastTriangles(std::vector<Triangle> & triangles,
                                 const int root,
                                 MPI_Comm comm)
  {
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    unsigned long long count =
        rank == root ? static_cast<unsigned long long>(triangles.size()) : 0;
    MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG_LONG, root, comm);

    nekrsCheck(count > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
                   count > static_cast<unsigned long long>(INT_MAX / 9),
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM triangle array cannot be represented safely for MPI broadcast\n");

    const std::size_t triangleCount = static_cast<std::size_t>(count);
    std::vector<double> packed(9 * triangleCount);
    if (rank == root)
      for (std::size_t triangle = 0; triangle < triangleCount; ++triangle)
        for (int vertex = 0; vertex < 3; ++vertex)
        {
          const std::size_t base = 9 * triangle + 3 * vertex;
          packed[base + 0] = triangles[triangle].vertex[vertex].x;
          packed[base + 1] = triangles[triangle].vertex[vertex].y;
          packed[base + 2] = triangles[triangle].vertex[vertex].z;
        }

    MPI_Bcast(packed.data(), static_cast<int>(packed.size()), MPI_DOUBLE, root, comm);

    if (rank != root)
    {
      triangles.resize(triangleCount);
      for (std::size_t triangle = 0; triangle < triangleCount; ++triangle)
      {
        triangles[triangle].sourceTriangle = triangle;
        for (int vertex = 0; vertex < 3; ++vertex)
        {
          const std::size_t base = 9 * triangle + 3 * vertex;
          triangles[triangle].vertex[vertex] = {
              packed[base + 0], packed[base + 1], packed[base + 2]};
        }
      }
    }
  }

  static std::uint64_t fnv1aAppend(std::uint64_t hash,
                                   const void * data,
                                   const std::size_t bytes)
  {
    const auto * values = static_cast<const unsigned char *>(data);
    for (std::size_t index = 0; index < bytes; ++index)
    {
      hash ^= values[index];
      hash *= UINT64_C(1099511628211);
    }
    return hash;
  }

  template <typename T>
  static std::uint64_t fnv1aVector(std::uint64_t hash, const std::vector<T> & values)
  {
    return fnv1aAppend(hash, values.data(), values.size() * sizeof(T));
  }

  static std::uint64_t markerHash(const MarkerSet & values)
  {
    std::uint64_t hash = UINT64_C(14695981039346656037);
    const std::uint64_t count = values.size();
    hash = fnv1aAppend(hash, &count, sizeof(count));
    hash = fnv1aVector(hash, values.coordinates);
    hash = fnv1aVector(hash, values.targetVelocity);
    hash = fnv1aVector(hash, values.volume);
    return fnv1aVector(hash, values.sourceTriangle);
  }

  void validateBroadcastResult(const MarkerSet & values, MPI_Comm comm)
  {
    const std::size_t count = values.size();
    bool localValid = count > 0 && values.coordinates.size() == 3 * count &&
                      values.targetVelocity.size() == 3 * count &&
                      values.sourceTriangle.size() == count;
    for (std::size_t p = 0; p < count && localValid; ++p)
    {
      localValid &= std::isfinite(values.volume[p]) && values.volume[p] > 0;
      localValid &= values.sourceTriangle[p] >= 0;
      for (int component = 0; component < 3; ++component)
      {
        localValid &= std::isfinite(values.coordinates[p + component * count]);
        localValid &= std::isfinite(values.targetVelocity[p + component * count]);
      }
    }

    int valid = localValid ? 1 : 0;
    int allValid = 0;
    MPI_Allreduce(&valid, &allValid, 1, MPI_INT, MPI_MIN, comm);
    nekrsCheck(!allValid,
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM received invalid marker data on at least one rank\n");

    const unsigned long long localHash = markerHash(values);
    unsigned long long minimumHash = 0;
    unsigned long long maximumHash = 0;
    MPI_Allreduce(&localHash, &minimumHash, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, comm);
    MPI_Allreduce(&localHash, &maximumHash, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm);
    nekrsCheck(minimumHash != maximumHash,
               comm,
               EXIT_FAILURE,
               "%s",
               "IBM marker arrays differ across MPI ranks\n");
    _broadcastValidationPassed = true;
  }

  void printGeometrySummary(const std::vector<Triangle> & triangles,
                            const SamplingStatistics & statistics) const
  {
    const double totalArea = std::accumulate(
        triangles.begin(), triangles.end(), 0.0,
        [](const double sum, const Triangle & triangle) {
          return sum + triangleArea(triangle);
        });
    const auto volumeRange = std::minmax_element(_markers.volume.begin(), _markers.volume.end());
    const double markerVolume =
        std::accumulate(_markers.volume.begin(), _markers.volume.end(), 0.0);
    printf("IBM geometry: triangles=%zu markers=%zu ranks=%d ell_marker=%.8e "
           "shellThickness=%.8e area=%.16e markerVolume=%.16e "
           "dVp=[%.8e, %.8e] leafMaxEdge=[%.8e, %.8e] "
           "stlTranslation=(%.8e, %.8e, %.8e)\n",
           triangles.size(),
           _markers.size(),
           platform->comm.mpiCommSize(),
           _markerSpacing,
           _shellThickness,
           totalArea,
           markerVolume,
           static_cast<double>(*volumeRange.first),
           static_cast<double>(*volumeRange.second),
           statistics.minimumLeafMaxEdge,
           statistics.maximumLeafMaxEdge,
           _stlTranslation[0],
           _stlTranslation[1],
           _stlTranslation[2]);
  }

  void buildDeviceData(const std::vector<Triangle> & triangles, MPI_Comm comm)
  {
    auto mesh = _nrs->meshV;
    std::vector<dfloat> x(mesh->Nlocal);
    std::vector<dfloat> y(mesh->Nlocal);
    std::vector<dfloat> z(mesh->Nlocal);
    std::vector<dfloat> jw(mesh->Nlocal);
    if (mesh->Nlocal > 0)
    {
      mesh->o_x.copyTo(x.data(), mesh->Nlocal);
      mesh->o_y.copyTo(y.data(), mesh->Nlocal);
      mesh->o_z.copyTo(z.data(), mesh->Nlocal);
      mesh->o_Jw.copyTo(jw.data(), mesh->Nlocal);
    }

    _maps = buildInteractionMaps(
        _markers, x, y, z, jw, _supportRadius, _gaussianWidth, comm, _periodicDimensions);

    _oInterpOffsets = deviceMemory<dlong>(_maps.interp.offsets);
    _oInterpIndices = deviceMemory<dlong>(_maps.interp.indices);
    _oInterpWeights = deviceMemory<dfloat>(_maps.interp.weights);
    _oSpreadOffsets = deviceMemory<dlong>(_maps.spread.offsets);
    _oSpreadIndices = deviceMemory<dlong>(_maps.spread.indices);
    _oSpreadWeights = deviceMemory<dfloat>(_maps.spread.weights);

    _oTargetVelocity = deviceMemory<dfloat>(_markers.targetVelocity);
    _oMarkerAcceleration = deviceMemory<dfloat>(_markers.targetVelocity);
    _oMarkerVelocity = deviceMemory<dfloat>(_markers.targetVelocity);
    _hostMarkerVelocity.resize(_markers.targetVelocity.size());

    if (_classifyElements)
      buildElementClassification(triangles, x, y, z, comm);
  }

  void buildElementClassification(const std::vector<Triangle> & triangles,
                                  const std::vector<dfloat> & x,
                                  const std::vector<dfloat> & y,
                                  const std::vector<dfloat> & z,
                                  MPI_Comm comm)
  {
    auto mesh = _nrs->meshV;
    const auto localRegions = classifyElements(triangles,
                                               x,
                                               y,
                                               z,
                                               static_cast<std::size_t>(mesh->Nelements),
                                               static_cast<std::size_t>(mesh->Np));
    const auto localCounts = countElementRegions(localRegions);
    const unsigned long long send[3] = {
        localCounts.fluid, localCounts.cut, localCounts.solid};
    unsigned long long recv[3] = {0, 0, 0};
    MPI_Allreduce(send, recv, 3, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    _globalRegionCounts.fluid = recv[0];
    _globalRegionCounts.cut = recv[1];
    _globalRegionCounts.solid = recv[2];

    const auto hostField =
        makeElementRegionField(localRegions, static_cast<std::size_t>(mesh->Np));
    _oElementRegionField = deviceMemory<dfloat>(hostField);

    if (platform->comm.mpiRank() == 0)
      printf("IBM element regions: fluid=%llu cut=%llu solid=%llu "
             "scalar00=(Fluid=0, Cut=1, Solid=2)\n",
             _globalRegionCounts.fluid,
             _globalRegionCounts.cut,
             _globalRegionCounts.solid);
  }

  void interpolate(const occa::memory & velocity)
  {
    ibm_interpolate(markerCount(),
                    _nrs->fluid->fieldOffset,
                    _oInterpOffsets,
                    _oInterpIndices,
                    _oInterpWeights,
                    velocity,
                    _oMarkerVelocity);
    _oMarkerVelocity.copyTo(_hostMarkerVelocity);
    MPI_Allreduce(MPI_IN_PLACE,
                  _hostMarkerVelocity.data(),
                  static_cast<int>(_hostMarkerVelocity.size()),
                  MPI_DFLOAT,
                  MPI_SUM,
                  platform->comm.mpiComm());
    _oMarkerVelocity.copyFrom(_hostMarkerVelocity);
  }

  nrs_t * _nrs = nullptr;
  MarkerSet _markers;
  InteractionMaps _maps;

  std::string _stlFile = "geometry/reference.stl";
  std::array<double, 3> _stlTranslation = {{0.0, 0.0, 0.0}};
  double _markerSpacing = 0.4;
  double _shellThickness = 0.05;
  double _supportRadius = 0.5;
  double _gaussianWidth = 0.25;
  int _maxSubdivisionDepth = 30;
  std::string _periodicDimensions = "none";
  bool _classifyElements = false;
  int _ciMode = 0;
  bool _broadcastValidationPassed = false;
  ElementRegionCounts _globalRegionCounts;

  deviceMemory<dlong> _oInterpOffsets;
  deviceMemory<dlong> _oInterpIndices;
  deviceMemory<dfloat> _oInterpWeights;
  deviceMemory<dlong> _oSpreadOffsets;
  deviceMemory<dlong> _oSpreadIndices;
  deviceMemory<dfloat> _oSpreadWeights;
  deviceMemory<dfloat> _oTargetVelocity;
  deviceMemory<dfloat> _oMarkerAcceleration;
  deviceMemory<dfloat> _oMarkerVelocity;
  deviceMemory<dfloat> _oElementRegionField;
  std::vector<dfloat> _hostMarkerVelocity;
};

} // namespace ibm
