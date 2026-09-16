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
