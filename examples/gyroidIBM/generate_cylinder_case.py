#!/usr/bin/env python3
"""Generate a structured cylinder O-grid mesh for nekRS enclosing the gyroid STL."""

from math import cos, pi, sin
from pathlib import Path
import struct
import numpy as np


ROOT = Path(__file__).resolve().parent
GEOMETRY = ROOT / "geometry"


def generate_cylinder_re2(
    path,
    xc=0.675,
    yc=0.675,
    R=0.561,
    z_min=0.0,
    z_max=2.70,
    nc=6,
    nr=4,
    nz=24,
):
    """
    Generate an O-grid (butterfly) hexahedral mesh for a circular cylinder.
    
    Parameters:
      path: output .re2 file path
      xc, yc: center of the cylinder axis
      R: outer cylinder radius (matches gyroid radius ~0.5553)
      z_min, z_max: streamwise domain extent
      nc: number of elements along each edge of the central square block
      nr: number of radial elements in each of the 4 outer quadrant blocks
      nz: number of elements in the streamwise z direction
    """
    # Half-width of central square block
    a = R * 0.45
    x_c = np.linspace(-a, a, nc + 1)
    y_c = np.linspace(-a, a, nc + 1)

    quads_2d = []
    outer_quad_indices = set()

    # 1. Central square block (nc x nc elements)
    for j in range(nc):
        for i in range(nc):
            p0 = (x_c[i], y_c[j])
            p1 = (x_c[i + 1], y_c[j])
            p2 = (x_c[i + 1], y_c[j + 1])
            p3 = (x_c[i], y_c[j + 1])
            quads_2d.append((p0, p1, p2, p3))

    # 2. 4 Outer quadrant blocks (East, North, West, South)
    def make_outer_block(inner_p0, inner_p1, theta0, theta1):
        grid = np.zeros((nr + 1, nc + 1, 2))
        for j in range(nc + 1):
            eta = j / nc
            theta = theta0 + (theta1 - theta0) * eta
            p_in = inner_p0 * (1 - eta) + inner_p1 * eta
            p_out = np.array([R * cos(theta), R * sin(theta)])
            for ir in range(nr + 1):
                xi = ir / nr
                grid[ir, j] = p_in * (1 - xi) + p_out * xi
        return grid

    g_east  = make_outer_block(np.array([a, -a]), np.array([a, a]), -pi / 4, pi / 4)
    g_north = make_outer_block(np.array([a, a]), np.array([-a, a]), pi / 4, 3 * pi / 4)
    g_west  = make_outer_block(np.array([-a, a]), np.array([-a, -a]), 3 * pi / 4, 5 * pi / 4)
    g_south = make_outer_block(np.array([-a, -a]), np.array([a, -a]), 5 * pi / 4, 7 * pi / 4)

    for g in [g_east, g_north, g_west, g_south]:
        for ir in range(nr):
            for j in range(nc):
                p0 = tuple(g[ir, j])
                p1 = tuple(g[ir + 1, j])
                p2 = tuple(g[ir + 1, j + 1])
                p3 = tuple(g[ir, j + 1])
                idx = len(quads_2d)
                quads_2d.append((p0, p1, p2, p3))
                if ir == nr - 1:
                    outer_quad_indices.add(idx)

    n_quads = len(quads_2d)
    element_count = n_quads * nz
    z_nodes = np.linspace(z_min, z_max, nz + 1)

    header = f"#v003{element_count:9d}{3:3d}{element_count:9d} gyroid cylinder"

    # Boundary ID convention for gyroidIBM.par:
    # boundaryTypeMap = udfDirichlet, zeroNeumann
    # ID 1: udfDirichlet (inlet at z=0 and no-slip wall at r=R)
    # ID 2: zeroNeumann (outflow at z=z_max)
    DIRICHLET_ID = 1
    NEUMANN_ID = 2
    boundaries = []

    with open(path, "wb") as stream:
        stream.write(header.encode("ascii").ljust(80, b" "))
        stream.write(struct.pack("<f", 6.54321))

        elem_id = 0
        for k in range(nz):
            z0 = z_nodes[k]
            z1 = z_nodes[k + 1]
            for q_idx, (p0, p1, p2, p3) in enumerate(quads_2d):
                elem_id += 1

                # 8 vertices of hexahedron (Nek5000 ordering)
                x = (p0[0] + xc, p1[0] + xc, p2[0] + xc, p3[0] + xc,
                     p0[0] + xc, p1[0] + xc, p2[0] + xc, p3[0] + xc)
                y = (p0[1] + yc, p1[1] + yc, p2[1] + yc, p3[1] + yc,
                     p0[1] + yc, p1[1] + yc, p2[1] + yc, p3[1] + yc)
                z = (z0, z0, z0, z0, z1, z1, z1, z1)

                stream.write(struct.pack("<d", 0.0))
                stream.write(struct.pack("<24d", *(x + y + z)))

                # Boundary conditions:
                # Inlet at z = z_min (face 5 in Nek convention: z- bottom)
                if k == 0:
                    boundaries.append((elem_id, 5, "v", DIRICHLET_ID))
                # Outlet at z = z_max (face 6 in Nek convention: z+ top)
                if k == nz - 1:
                    boundaries.append((elem_id, 6, "O", NEUMANN_ID))
                # Outer cylinder wall at r = R (face 2 in Nek convention: +x local)
                if q_idx in outer_quad_indices:
                    boundaries.append((elem_id, 2, "v", DIRICHLET_ID))

        # Curved sides (straight elements: ncurve = 0)
        stream.write(struct.pack("<d", 0.0))

        # Boundary conditions
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

    print(f"Successfully generated {path}:")
    print(f"  Elements: {element_count} ({n_quads} quads/layer x {nz} layers)")
    print(f"  Center: ({xc}, {yc}), Radius: {R}")
    print(f"  Streamwise extent: [{z_min}, {z_max}]")
    print(f"  Boundaries: {len(boundaries)} (inlet=5, outlet=6, wall=2)")


if __name__ == "__main__":
    out_file = ROOT / "gyroidIBM_cylinder.re2"
    generate_cylinder_re2(out_file)
