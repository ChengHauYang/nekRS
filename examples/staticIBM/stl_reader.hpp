#pragma once

#include "ibmGeometry.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ibm
{
namespace detail
{

inline std::uint32_t decodeLittleEndianU32(const unsigned char * bytes)
{
  return static_cast<std::uint32_t>(bytes[0]) |
         (static_cast<std::uint32_t>(bytes[1]) << 8) |
         (static_cast<std::uint32_t>(bytes[2]) << 16) |
         (static_cast<std::uint32_t>(bytes[3]) << 24);
}

inline float decodeLittleEndianFloat(const unsigned char * bytes)
{
  const std::uint32_t bits = decodeLittleEndianU32(bytes);
  float value;
  static_assert(sizeof(value) == sizeof(bits), "binary STL requires 32-bit float");
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

inline void readExact(std::ifstream & stream,
                      char * destination,
                      const std::streamsize bytes,
                      const std::string & context)
{
  stream.read(destination, bytes);
  if (stream.gcount() != bytes || !stream)
    throw std::runtime_error("failed to read " + context + " from binary STL");
}

} // namespace detail

inline std::vector<Triangle> readBinaryStl(const std::string & path)
{
  std::ifstream stream(path, std::ios::binary | std::ios::ate);
  if (!stream)
    throw std::runtime_error("cannot open binary STL file: " + path);

  const std::streamoff fileSize = stream.tellg();
  if (fileSize < 84)
    throw std::runtime_error("binary STL is shorter than its 84-byte header: " + path);
  stream.seekg(0, std::ios::beg);

  std::array<unsigned char, 84> header{};
  detail::readExact(stream,
                    reinterpret_cast<char *>(header.data()),
                    static_cast<std::streamsize>(header.size()),
                    "header");
  const std::uint32_t triangleCount = detail::decodeLittleEndianU32(header.data() + 80);
  if (triangleCount == 0)
    throw std::runtime_error("binary STL contains zero triangles: " + path);

  constexpr std::uint64_t headerBytes = 84;
  constexpr std::uint64_t recordBytes = 50;
  if (triangleCount > (std::numeric_limits<std::uint64_t>::max() - headerBytes) / recordBytes)
    throw std::overflow_error("binary STL triangle count causes a size overflow");
  const std::uint64_t expectedSize = headerBytes + recordBytes * triangleCount;
  if (fileSize < 0 || static_cast<std::uint64_t>(fileSize) != expectedSize)
  {
    std::ostringstream message;
    message << "binary STL size mismatch for " << path << ": expected " << expectedSize
            << " bytes, found " << fileSize;
    throw std::runtime_error(message.str());
  }

  std::vector<Triangle> triangles;
  triangles.reserve(triangleCount);
  Vec3 lower{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
  Vec3 upper{-std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};

  std::array<unsigned char, recordBytes> record{};
  for (std::uint32_t index = 0; index < triangleCount; ++index)
  {
    detail::readExact(stream,
                      reinterpret_cast<char *>(record.data()),
                      static_cast<std::streamsize>(record.size()),
                      "triangle record " + std::to_string(index));
    Triangle triangle{};
    triangle.sourceTriangle = index;
    // Bytes 0..11 store an untrusted normal. Vertices start at byte 12.
    for (unsigned int vertex = 0; vertex < 3; ++vertex)
    {
      const unsigned char * data = record.data() + 12 + 12 * vertex;
      triangle.vertex[vertex] = {
          static_cast<double>(detail::decodeLittleEndianFloat(data)),
          static_cast<double>(detail::decodeLittleEndianFloat(data + 4)),
          static_cast<double>(detail::decodeLittleEndianFloat(data + 8))};
      const Vec3 & value = triangle.vertex[vertex];
      if (!detail::finite(value))
        throw std::runtime_error("binary STL triangle " + std::to_string(index) +
                                 " contains a non-finite vertex");
      lower.x = std::min(lower.x, value.x);
      lower.y = std::min(lower.y, value.y);
      lower.z = std::min(lower.z, value.z);
      upper.x = std::max(upper.x, value.x);
      upper.y = std::max(upper.y, value.y);
      upper.z = std::max(upper.z, value.z);
    }
    triangles.push_back(triangle);
  }

  const double diagonal = detail::edgeLength(lower, upper);
  const double areaTolerance =
      256.0 * std::numeric_limits<double>::epsilon() * diagonal * diagonal;
  for (std::size_t index = 0; index < triangles.size(); ++index)
  {
    const double area = detail::triangleArea(triangles[index]);
    if (!std::isfinite(area) || area <= areaTolerance)
    {
      std::ostringstream message;
      message << "binary STL triangle " << index << " is degenerate: area=" << area
              << ", scale-aware tolerance=" << areaTolerance;
      throw std::runtime_error(message.str());
    }
  }
  return triangles;
}

} // namespace ibm
