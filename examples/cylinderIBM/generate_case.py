#!/usr/bin/env python3
"""Generate the deterministic cylinder STL and one-element-span Nek mesh."""

from math import cos, pi, sin
from pathlib import Path
import struct


ROOT = Path(__file__).resolve().parent
GEOMETRY = ROOT / "geometry"


def write_triangle(stream, vertices):
    a, b, c = vertices
    ab = tuple(b[i] - a[i] for i in range(3))
    ac = tuple(c[i] - a[i] for i in range(3))
    normal = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    magnitude = sum(value * value for value in normal) ** 0.5
    stream.write(struct.pack("<3f", *(value / magnitude for value in normal)))
    for vertex in vertices:
        stream.write(struct.pack("<3f", *vertex))
    stream.write(struct.pack("<H", 0))


def write_cylinder(path, radius=0.5, length=1.0, n_theta=64, n_z=8):
    triangles = []
    for k in range(n_z):
        z0 = length * k / n_z
        z1 = length * (k + 1) / n_z
        for sector in range(n_theta):
            theta0 = 2.0 * pi * sector / n_theta
            theta1 = 2.0 * pi * (sector + 1) / n_theta
            v00 = (radius * cos(theta0), radius * sin(theta0), z0)
            v10 = (radius * cos(theta1), radius * sin(theta1), z0)
            v01 = (radius * cos(theta0), radius * sin(theta0), z1)
            v11 = (radius * cos(theta1), radius * sin(theta1), z1)
            triangles.extend(((v00, v10, v01), (v10, v11, v01)))

    header = b"nekRS cylinderIBM radius=0.5 z-periodic length=1.0 no end caps"
    with path.open("wb") as stream:
        stream.write(header.ljust(80, b"\0"))
        stream.write(struct.pack("<I", len(triangles)))
        for triangle in triangles:
            write_triangle(stream, triangle)


def write_re2(path, nx=48, ny=32, nz=1):
    x_min, x_max = -8.0, 16.0
    y_min, y_max = -8.0, 8.0
    z_min, z_max = 0.0, 1.0
    element_count = nx * ny * nz
    header = f"#v003{element_count:9d}{3:3d}{element_count:9d} cylinderIBM"

    with path.open("wb") as stream:
        stream.write(header.encode("ascii").ljust(80, b" "))
        stream.write(struct.pack("<f", 6.54321))

        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    x0 = x_min + (x_max - x_min) * i / nx
                    x1 = x_min + (x_max - x_min) * (i + 1) / nx
                    y0 = y_min + (y_max - y_min) * j / ny
                    y1 = y_min + (y_max - y_min) * (j + 1) / ny
                    z0 = z_min + (z_max - z_min) * k / nz
                    z1 = z_min + (z_max - z_min) * (k + 1) / nz
                    x = (x0, x1, x1, x0, x0, x1, x1, x0)
                    y = (y0, y0, y1, y1, y0, y0, y1, y1)
                    z = (z0, z0, z0, z0, z1, z1, z1, z1)
                    stream.write(struct.pack("<d", 0.0))
                    stream.write(struct.pack("<24d", *(x + y + z)))

        stream.write(struct.pack("<d", 0.0))

        # boundaryTypeMap in cylinderIBM.par is (udfDirichlet, zeroNeumann),
        # so boundary IDs are: 1 = udfDirichlet (freestream inlet + top/bottom),
        # 2 = zeroNeumann (outflow). Periodic faces carry no ID.
        DIRICHLET_ID = 1
        NEUMANN_ID = 2
        boundaries = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    element = 1 + i + nx * (j + ny * k)
                    if j == 0:
                        boundaries.append((element, 1, "v", 0, 0, DIRICHLET_ID))
                    if i == nx - 1:
                        boundaries.append((element, 2, "O", 0, 0, NEUMANN_ID))
                    if j == ny - 1:
                        boundaries.append((element, 3, "v", 0, 0, DIRICHLET_ID))
                    if i == 0:
                        boundaries.append((element, 4, "v", 0, 0, DIRICHLET_ID))
                    if k == 0:
                        peer = element + nx * ny * (nz - 1)
                        boundaries.append((element, 5, "P", peer, 6, 0))
                    if k == nz - 1:
                        peer = element - nx * ny * (nz - 1)
                        boundaries.append((element, 6, "P", peer, 5, 0))

        stream.write(struct.pack("<d", float(len(boundaries))))
        for element, side, code, peer_element, peer_side, bid in boundaries:
            # bc(5) — the 5th BC parameter (7th value in the packed record) —
            # carries the integer boundary ID that nekRS reads via
            # nekInterface.f (boundaryID = bc(5,ifc,iel,ifld_bId)). Leaving it
            # zero causes nekRS to report NboundaryIDs=0 and abort.
            stream.write(
                struct.pack(
                    "<7d",
                    float(element),
                    float(side),
                    float(peer_element),
                    float(peer_side),
                    0.0,
                    0.0,
                    float(bid),
                )
            )
            stream.write(f"{code:<3}".encode("ascii").ljust(8, b" "))


if __name__ == "__main__":
    GEOMETRY.mkdir(exist_ok=True)
    write_cylinder(GEOMETRY / "cylinder.stl")
    write_re2(ROOT / "cylinderIBM.re2")
    print("generated cylinderIBM.re2 and geometry/cylinder.stl")
