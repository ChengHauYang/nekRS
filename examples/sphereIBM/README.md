# Lid-driven cavity with an immersed sphere

This case exercises the static direct-forcing IBM in a three-dimensional
lid-driven cavity. The configuration is

- domain: `[0, 1]^3`;
- stationary sphere: center `(0.5, 0.5, 0.5)`, radius `0.25`;
- moving top wall (`y = 1`): velocity `(1, 0, 0)`;
- other cavity walls: stationary no-slip;
- density: `1`;
- kinematic viscosity: `0.01`.

Using the cavity side length and lid velocity, the Reynolds number is
`Re = U L / nu = 100`.

Run the case with

```console
mpirun -np 2 nekrs --setup sphereIBM.par
```

The case reports the maximum and RMS marker slip and the integrated IBM force
at checkpoints and at the final step. `generate_case.py` deterministically
regenerates the binary sphere STL and the structured `8 x 8 x 8` Nek mesh.
