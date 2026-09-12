#!/usr/bin/env python3
"""Generate the deterministic sphere STL and structured unit-cube Nek mesh."""

from math import cos, pi, sin
from pathlib import Path
import struct


ROOT = Path(__file__).resolve().parent
GEOMETRY = ROOT / "geometry"


def write_sphere(path, center=(0.5, 0.5, 0.5), radius=0.25, n_latitude=20, n_longitude=40):
    vertices = [(center[0], center[1] + radius, center[2])]
    for latitude in range(1, n_latitude):
        phi = pi * latitude / n_latitude
        for longitude in range(n_longitude):
            theta = 2 * pi * longitude / n_longitude
            vertices.append(
                (
                    center[0] + radius * sin(phi) * cos(theta),
                    center[1] + radius * cos(phi),
                    center[2] + radius * sin(phi) * sin(theta),
                )
            )
    south = len(vertices)
    vertices.append((center[0], center[1] - radius, center[2]))

    faces = []
    first_ring = 1
    for longitude in range(n_longitude):
        next_longitude = (longitude + 1) % n_longitude
        faces.append((0, first_ring + next_longitude, first_ring + longitude))

    for latitude in range(n_latitude - 2):
        upper = 1 + latitude * n_longitude
        lower = upper + n_longitude
        for longitude in range(n_longitude):
            next_longitude = (longitude + 1) % n_longitude
            faces.append((upper + longitude, upper + next_longitude, lower + longitude))
            faces.append((upper + next_longitude, lower + next_longitude, lower + longitude))

    last_ring = 1 + (n_latitude - 2) * n_longitude
    for longitude in range(n_longitude):
        next_longitude = (longitude + 1) % n_longitude
        faces.append((last_ring + longitude, last_ring + next_longitude, south))

    header = b"nekRS sphereIBM center=(0.5,0.5,0.5) radius=0.25"
    with path.open("wb") as stream:
        stream.write(header.ljust(80, b"\0"))
        stream.write(struct.pack("<I", len(faces)))
        for face in faces:
            a, b, c = (vertices[index] for index in face)
            ab = tuple(b[i] - a[i] for i in range(3))
            ac = tuple(c[i] - a[i] for i in range(3))
            normal = (
                ab[1] * ac[2] - ab[2] * ac[1],
                ab[2] * ac[0] - ab[0] * ac[2],
                ab[0] * ac[1] - ab[1] * ac[0],
            )
            magnitude = sum(value * value for value in normal) ** 0.5
            stream.write(struct.pack("<3f", *(value / magnitude for value in normal)))
            for vertex in (a, b, c):
                stream.write(struct.pack("<3f", *vertex))
            stream.write(struct.pack("<H", 0))


def write_unit_cube_re2(path, elements_per_direction=8):
    n = elements_per_direction
    element_count = n**3
    header = f"#v003{element_count:9d}{3:3d}{element_count:9d} sphereIBM unit cube"

    with path.open("wb") as stream:
        stream.write(header.encode("ascii").ljust(80, b" "))
        stream.write(struct.pack("<f", 6.54321))

        for k in range(n):
            for j in range(n):
                for i in range(n):
                    x0, x1 = i / n, (i + 1) / n
                    y0, y1 = j / n, (j + 1) / n
                    z0, z1 = k / n, (k + 1) / n
                    x = (x0, x1, x1, x0, x0, x1, x1, x0)
                    y = (y0, y0, y1, y1, y0, y0, y1, y1)
                    z = (z0, z0, z0, z0, z1, z1, z1, z1)
                    stream.write(struct.pack("<d", 0.0))
                    stream.write(struct.pack("<24d", *(x + y + z)))

        # No curved sides.
        stream.write(struct.pack("<d", 0.0))

        boundaries = []
        for k in range(n):
            for j in range(n):
                for i in range(n):
                    element = 1 + i + n * (j + n * k)
                    if j == 0:
                        boundaries.append((element, 1))
                    if i == n - 1:
                        boundaries.append((element, 2))
                    if j == n - 1:
                        boundaries.append((element, 3))
                    if i == 0:
                        boundaries.append((element, 4))
                    if k == 0:
                        boundaries.append((element, 5))
                    if k == n - 1:
                        boundaries.append((element, 6))

        stream.write(struct.pack("<d", float(len(boundaries))))
        for element, side in boundaries:
            # bc(5) carries the integer boundary ID that nekRS reads via
            # nekInterface.f (boundaryID = bc(5,ifc,iel,ifld_bId)); leaving it
            # zero causes nekRS to see NboundaryIDs=0. All six outer faces
            # share one ID because sphereIBM.par declares a single
            # boundaryTypeMap entry (udfDirichlet), and the udf discriminates
            # top-wall vs. no-slip by y-coordinate rather than by ID.
            stream.write(struct.pack("<7d", float(element), float(side), 0, 0, 0, 0, 1.0))
            stream.write(b"W  ".ljust(8, b" "))


if __name__ == "__main__":
    GEOMETRY.mkdir(exist_ok=True)
    write_sphere(GEOMETRY / "sphere_r025.stl")
    write_unit_cube_re2(ROOT / "sphereIBM.re2")
    print("generated sphereIBM.re2 and geometry/sphere_r025.stl")
