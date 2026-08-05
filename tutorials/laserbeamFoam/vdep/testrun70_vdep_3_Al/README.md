# testrun70 — beam quality remedy (Zixun)

This case is **config-only** — `constant/polyMesh` and a seeded timestep are deliberately not
included (that's large binary data, better generated locally than pushed through git).

Already done here:
- All of testrun64's config copied over (`constant/`, `system/`) — same 750W, endTime=400µs,
  decomposition, etc.
- The one remedy change applied: `constant/LaserProperties`'s `Radius_Flavour` set to `0.7`
  (was `1.336`).
- Mass-conservation monitor added to `system/controlDict`'s `functions` block.

**Still needed before running** — the seed hand-off from testrun61 (see
[vdep_remedy_sims.md](../../../../vdep_remedy_sims.md), step 3 of "Per-case setup"):
1. Reconstruct testrun61's latest timestep (if not already done).
2. Copy `testrun61_vdep_3_Al/constant/polyMesh` and its latest timestep directory into this case.
3. `decomposePar -time <TIME>` here (uses this case's own `decomposeParDict`, already copied).
4. `checkMesh -time <TIME>` — expect `Mesh OK` before running.

No solver changes needed for this case — runs on the plain `lbf3` image.
