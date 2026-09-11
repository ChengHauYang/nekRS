# Quasi-two-dimensional flow past an immersed cylinder

This case exercises periodic IBM interactions using a nominally two-dimensional
flow-past-a-cylinder benchmark. The solver remains three-dimensional, but the
mesh has one spectral element in the spanwise direction and periodic boundary
conditions on its two `z` faces.

The configuration is

- domain: `[-8, 16] x [-8, 8] x [0, 1]`;
- mesh: `48 x 32 x 1` spectral elements;
- stationary cylinder: diameter `D = 1`, centered at `(0, 0)`, axis parallel to `z`;
- inlet and lateral far field: velocity `(1, 0, 0)`;
- outlet: zero-Neumann velocity;
- `z = 0` and `z = 1`: periodic;
- density: `1`;
- kinematic viscosity: `0.01`;
- Reynolds number: `Re_D = U_infinity D / nu = 100`.

The STL contains only the lateral cylinder surface. It intentionally has no end
caps because the cylinder crosses the periodic boundary. The IBM interaction
map uses the nearest periodic image of mesh nodes when a marker support sphere
crosses `z = 0` or `z = 1`.

Regenerate the deterministic mesh and STL with

```console
python3 generate_case.py
```

Run the case with

```console
mpirun -np 2 nekrs --setup cylinderIBM.par
```

The case reports marker slip, integrated IBM force, and lift and drag
coefficients at checkpoints. Because this is a one-element spanwise extrusion,
the initial condition, boundary data, and geometry must remain independent of
`z`, and the spanwise velocity must remain zero.
