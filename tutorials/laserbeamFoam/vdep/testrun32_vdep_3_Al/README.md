# testrun32_vdep_3_Al — No Results (Superseded by testrun33)

## Why there are no results here

This case was run and reconstructed, but its results were discarded. The laser scan track
starts only **150 µm** from the beginning (z = 0) boundary of the domain. This proximity
to the inlet boundary can introduce artefacts: the melt pool dynamics near the start of
the track may be influenced by boundary effects rather than reflecting the true material
response to the laser.

## What was changed in testrun33

`testrun33_vdep_3_Al` is identical to this case in every other respect (500 W, 1.5 m/s,
35 µm beam radius, AlSi10Mg) but moves the laser starting position to **500 µm** from the
beginning edge, providing adequate clearance from the boundary. The end position remains
at 100 µm from the far edge (z = 2400 µm), consistent with testrun32.

| Parameter        | testrun32     | testrun33     |
|------------------|---------------|---------------|
| Laser start (z)  | 150 µm        | **500 µm**    |
| Laser end (z)    | 2400 µm       | 2400 µm       |
| Track length     | 2250 µm       | 1900 µm       |
| Simulation time  | 1.500 ms      | 1.267 ms      |

Use `testrun33_vdep_3_Al` for all analysis.
