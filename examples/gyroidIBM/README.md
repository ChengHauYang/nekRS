# gyroidIBM

This example runs flow through a gyroid STL immersed in a structured box mesh.
The solid is represented by the symmetric static-surface IBM implementation.

The box mesh spans `[0, 1.35] x [0, 1.35] x [0, 2.70]`. The STL is stored
losslessly as split `geometry/parts/gyroid.stl.xz.part*` files;
`generate_case.py` restores `geometry/gyroid.stl` for runs. The STL is kept in
its source coordinates and placed into the box through the input file:

```ini
[CASEDATA]
stl_file = geometry/gyroid.stl
stl_translation = 0.675, 0.675, 0.7575
classify_elements = true
```

`classify_elements = true` writes the element region marker as `scalar00` in
the checkpoint fields. Values are:

| Value | Region |
| ---: | --- |
| 0 | Fluid |
| 1 | Cut |
| 2 | Solid |

The inlet is the `z = 0` face, the outlet is the `z = 2.70` face, and the
remaining box faces are no-slip walls. The STL and mesh can be regenerated with:

```bash
python3 generate_case.py
```

The default mesh is intentionally modest for setup and inspection. Increase the
element counts in `generate_case.py` and tighten the IBM length scales in
`gyroidIBM.par` for production-quality runs.

## Build

This example uses a development-only IBM API (`hasElementClassification`,
`elementRegionField`, `elementRegionCounts`), so a build from the current
source tree (not a release install) is required:

```bash
cd ~/packages/nekRS
cmake --build build -j$(nproc)     # build (requires build/ configured with CUDA)
cmake --install build              # install (overwrites ~/.local/nekrs)
```

## Run

Requires the `nekrs` binary on your `PATH`. One MPI rank per GPU is the default
(`device-id = LOCAL-RANK`); with `-np` larger than the GPU count, add
`--device-id 0` so all ranks share the first GPU.

A single `--setup` invocation compiles the UDF and runs the simulation:

```bash
mpirun -np 2 nekrs --setup gyroidIBM.par            # 1 rank per GPU
# or force all ranks onto GPU 0:
mpirun -np 6 nekrs --setup gyroidIBM.par --device-id 0
```

Run from this directory so `gyroidIBM.par` and the case files resolve.
