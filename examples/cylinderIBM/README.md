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
- Reynolds number: `Re_D = U_infinity D / nu = 100`;
- final nondimensional time: `150`;
- statistics window: `80 <= t <= 150`.

The TOMBO3 scheme treats convection explicitly, so its stable step remains
limited by the convective CFL even though the pressure and viscous solves are
implicit. The case requests a target CFL of `0.5`, starts from `dt = 0.01`, and
caps adaptive time steps at `0.01`:

```ini
dt = targetCFL=0.5 + max=1e-2 + initial=1e-2
```

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

The built-in checkpoint schedule is disabled. Starting at `t = 80`, the UDF
writes full flow fields every `0.5` time units and reports marker slip,
integrated IBM force, and lift and drag coefficients every `0.1` time units.
This provides about 12 field snapshots and about 60 force samples per expected
vortex-shedding period at `Re_D = 100` without storing unnecessary transient
fields. Because this is a one-element spanwise extrusion, the initial condition,
boundary data, and geometry must remain independent of `z`, and the spanwise
velocity must remain zero.
