#!/usr/bin/env python3
"""Generate the structured box Nek mesh and restore the gyroid STL."""

import lzma
from pathlib import Path
import struct


ROOT = Path(__file__).resolve().parent
GEOMETRY = ROOT / "geometry"


def restore_gyroid_stl():
    stl_path = GEOMETRY / "gyroid.stl"
    compressed_path = GEOMETRY / "gyroid.stl.xz"
    if stl_path.exists():
        return
    if compressed_path.exists():
        compressed = compressed_path.read_bytes()
    else:
        parts = sorted((GEOMETRY / "parts").glob("gyroid.stl.xz.part*"))
        if not parts:
            raise FileNotFoundError("missing gyroid.stl.xz or gyroid.stl.xz.part files")
        compressed = b"".join(part.read_bytes() for part in parts)
    stl_path.write_bytes(lzma.decompress(compressed))


def write_box_re2(path, nx=8, ny=8, nz=16):
    x_min, x_max = 0.0, 1.35
    y_min, y_max = 0.0, 1.35
    z_min, z_max = 0.0, 2.70
    element_count = nx * ny * nz
    header = f"#v003{element_count:9d}{3:3d}{element_count:9d} gyroidIBM box"

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

        # boundaryTypeMap in gyroidIBM.par is (udfDirichlet, zeroNeumann).
        # IDs: 1 = inlet + box walls, 2 = outflow.
        DIRICHLET_ID = 1
        NEUMANN_ID = 2
        boundaries = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    element = 1 + i + nx * (j + ny * k)
                    if j == 0:
                        boundaries.append((element, 1, "v", DIRICHLET_ID))
                    if i == nx - 1:
                        boundaries.append((element, 2, "v", DIRICHLET_ID))
                    if j == ny - 1:
                        boundaries.append((element, 3, "v", DIRICHLET_ID))
                    if i == 0:
                        boundaries.append((element, 4, "v", DIRICHLET_ID))
                    if k == 0:
                        boundaries.append((element, 5, "v", DIRICHLET_ID))
                    if k == nz - 1:
                        boundaries.append((element, 6, "O", NEUMANN_ID))

        stream.write(struct.pack("<d", float(len(boundaries))))
        for element, side, code, boundary_id in boundaries:
            stream.write(
                struct.pack(
                    "<7d",
                    float(element),
                    float(side),
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    float(boundary_id),
                )
            )
            stream.write(f"{code:<3}".encode("ascii").ljust(8, b" "))


if __name__ == "__main__":
    restore_gyroid_stl()
    write_box_re2(ROOT / "gyroidIBM.re2")
    print("generated gyroidIBM.re2 and restored geometry/gyroid.stl")
