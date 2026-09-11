#!/usr/bin/env python3
"""Generate the deterministic four-triangle binary STL used by Phase 1 CI."""

from pathlib import Path
import struct


vertices = (
    (-0.25, -0.25, -0.25),
    (0.25, -0.25, -0.25),
    (-0.25, 0.25, -0.25),
    (-0.25, -0.25, 0.25),
)

# Fixed source-triangle order. Stored normals are deliberately zero because the
# reader must reconstruct geometry from vertices rather than trust STL normals.
faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))

output = Path(__file__).with_name("reference.stl")
header = b"nekRS staticIBM deterministic tetrahedron fixture"
header = header.ljust(80, b"\0")

with output.open("wb") as stream:
    stream.write(header)
    stream.write(struct.pack("<I", len(faces)))
    for face in faces:
        stream.write(struct.pack("<3f", 0.0, 0.0, 0.0))
        for vertex in face:
            stream.write(struct.pack("<3f", *vertices[vertex]))
        stream.write(struct.pack("<H", 0))

print(output)
