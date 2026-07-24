# Jul 23 Tasks — pinned geometry/mesh + power-sweep case set (650–900 W)

## Context (abstracted from user's request)

Starting a new set of `laserbeamFoam` VDEP case files built around a laser-power
sweep. All new cases will share one geometry and mesh, which needs to be pinned
down first (**Task 1**, this document). Once pinned, two seed cases will be built
following the tr50/tr51 pattern (staged run: startup stage, then a settle/continuation
stage) to serve as the common starting checkpoint. From that shared seed, six
production cases will branch, each run at a different constant laser power:
**650, 700, 750, 800, 850, 900 W**.

A planned departure from tr50–54: those cases used `dynamicRefineFvMesh`, an
AMR scheme that refines wherever the `alpha.metal` interface field currently
sits. Since tr50, tr51 (and their continuations tr52/tr53/tr54) already show
where "the action" is — the melt pool / keyhole stays close to a known surface
height and close to the laser track's centerline — the new plan is to replace
the dynamic, field-triggered AMR with a **static, location-based refinement**:
a fixed refined region defined by y (height/depth relative to the metal
surface) and lateral x-distance from the track centerline, rather than mesh
adaptation driven by the live `alpha.metal` field.

## Task 1 — Pin down shared geometry & mesh

Values below are as currently set in `tr50`/`tr51`
(`tutorials/laserbeamFoam/vdep/testrun50_vdep_3_Al`,
`tutorials/laserbeamFoam/vdep/testrun51_vdep_3_Al`), which are the reference
cases for the new set.

### Domain

`convertToMeters 0.001` (all vertex coordinates below are in mm). **Revised
2026-07-23**: gas headspace and metal thickness both clipped down from
tr50/tr51's values to concentrate the mesh budget nearer the surface —
domain height cut from 0.7mm to 0.5mm:

| Axis | Range | Extent | Meaning |
|---|---|---|---|
| x | -0.32 → 0.32 mm | 0.64 mm | lateral, across the track |
| y | 0 → 0.5 mm (was 0 → 0.7 mm) | 0.5 mm (was 0.7 mm) | vertical (height/depth) |
| z | 0 → 3.2 mm | 3.2 mm | scan / build direction |

Base mesh: `hex (16 13 80) simpleGrading (1 1 1)` → cell size ≈ 40 µm (x) ×
38.5 µm (y, = 500/13) × 40 µm (z). Base cell count: 16×13×80 = **16,640**
(down from 23,040 — a ~28% cut from the y-clip alone, before any refinement
is even applied).

Boundary patches: `atmosphere` at y=0 (open top), `lowerWall` at y=0.5mm
(solid substrate, type `wall`), `rightWall`/`leftWall` at x=±0.32mm,
`frontAndBack` at z=0 and z=3.2mm.

### Where the metal surface is

`setFieldsDict` seeds `alpha.metal = 1` inside the box
`y ∈ [0.2mm, 0.5mm]` (everywhere else defaults to 0 = gas):

- Metal slab: y = 0.2 → 0.5 mm (**0.3 mm thick**, was 0.4mm, down to the substrate/`lowerWall`)
- Gas headspace: y = 0 → 0.2 mm (**0.2 mm**, was 0.3mm, up to the `atmosphere` boundary)
- **Metal free surface (melt-pool interface) sits at y = 0.2 mm** (was 0.3mm).

`V_incident (0 1 0)` — the laser travels in +y, i.e. from the `atmosphere`
boundary (y=0) down onto the metal surface (y=0.2mm).

### Laser track

From `constant/LaserProperties` and `constant/timeVsLaserPosition`, adjusted
for the new surface height:

- `laserRadius` = 3.5e-5 m (35 µm), `nRadial` 10, `nAngular` 36, `wavelength` 1.064e-6 m
- Track is a straight line along z, at fixed x = 0 (domain centerline) and
  fixed y = 0.2 mm (on the new metal surface):
  - t = 0 → (x=0, y=0.2mm, z=0.5mm)
  - t = 4e-4 s → (x=0, y=0.2mm, z=2.9mm)
- Track length: 2.9 − 0.5 = 2.4 mm over 400 µs → **scan speed 6 m/s**
- Margins: 0.5 mm clearance from the z=0 start face, 0.3 mm clearance to the
  z=3.2mm end face

### Maximum cells & current refinement method (dynamic AMR, to be replaced)

`constant/dynamicMeshDict` (tr50/tr51 values):

- `dynamicFvMesh dynamicRefineFvMesh`, `refineInterval 5`
- Refinement field: `alpha.metal`, refine where 0.001 < α < 0.999 (the
  metal/gas interface band), unrefine below 0.001
- `maxRefinement 3` → 3 doublings, finest cell ≈ 40µm / 8 = **5 µm**
- `nBufferLayers 3`, `maxCells 800,000` (base mesh is 23,040 cells; 3-level
  AMR near the interface adds roughly 150k–300k cells in practice)
- `dumpLevel true` (writes the refinement-level field for inspection)

**2026-07-23 feedback from user:** in prior runs the simulation hit the
`maxCells` ceiling and the last ~30% of the run was left under-refined
(coarse cells in regions that still needed refinement) — `maxCells=800,000`
was not enough. This is a hard constraint on the new case set: whatever
budget is set (static box + AMR headroom) must not repeat this, i.e. the
ceiling needs to be sized with margin for the whole 400µs run, not just the
early stage.

### Proposed change: static, location-based refinement

Instead of (or in addition to) triggering refinement dynamically off the
live `alpha.metal` field, define a fixed refined region up front, using:

- a **y-band** (height/depth) around the known metal surface / melt-pool
  extent, and
- an **x-band** (lateral distance from the track centerline, x=0).

**Candidate box (2026-07-23, relative to the new y=0.2mm surface):**
y ∈ [surface − 0.1mm, surface + 0.25mm] = [0.1mm, 0.45mm] (100µm above the
surface, 250µm below it, 350µm total — leaves 50µm of metal clearance above
the `lowerWall` at 0.5mm), x ∈ [−0.2mm, 0.2mm] (400µm wide), z = full domain
(3.2mm, unrestricted).

Base cells covered by this box: 10 (x) × 9 (y) × 80 (z) = **7,200** of the
16,640 total base cells.

| Target resolution | Levels | Subcells/base cell | Cells in box |
|---|---|---|---|
| 20 µm | 1 | 8 | 57,600 |
| 10 µm | 2 | 64 | 460,800 |
| 5 µm | 3 | 512 | 3,686,400 |

Full-z coverage at 5µm (matching tr50/51's finest AMR level) is ~3.7M cells
— 4.6–24× more than tr50/51's actual practical range (150k–300k) and well
past their `maxCells=800,000`, which the user has now confirmed was already
too low. 10µm (460,800 cells) is the more realistic static target.

**Recommendation to keep adaptive capability:** don't discard
`dynamicRefineFvMesh` — point it at a static geometric indicator field
(1 inside the box above, 0 outside; optionally combined with `alpha.metal`
via `max()`) instead of the raw interface field. This keeps the existing
grading/buffer-layer machinery, makes refinement deterministic in the known
action zone, and leaves `maxCells` headroom for AMR to still push to 5µm
wherever `alpha.metal` demands it outside the box — including late in the
run, which is exactly where tr50–54 ran out of budget before.

**Revised box (2026-07-23, tightened lateral band):** y ∈ [0.1mm, 0.4mm]
(100µm above surface, 200µm below — 300µm total), x ∈ [−0.1mm, 0.1mm]
(200µm wide, half the previous lateral band), z = full domain. Estimated
~3,120 base cells in box → ~1.61M cells at 5µm, ~213k at 10µm.

User chose **5µm** for this box (option B) and asked for confirmation that
`maxCells` still leaves `dynamicRefineFvMesh` headroom to refine outside the
box in any direction if the melt pool/keyhole/spatter goes beyond it.
Recommended `maxCells = 2,500,000` (~890k cells of headroom above the
static floor).

### Mesh-generation test (`testrun55_vdep_mesh_test`, 2026-07-23)

Built the case on the pinned geometry above (`blockMesh` → `hex (16 13 80)`
= 16,640 base cells) and applied the box refinement via 3 passes of
`topoSet` (boxToCell, `system/topoSetDict`) + `refineHexMesh -overwrite`
(mesh-generation only, **no solver run**). Actual results vs. the estimate:

| Stage | Cells | Estimate |
|---|---|---|
| Base mesh | 16,640 | 16,640 (exact) |
| After level 1 (20µm) | 37,640 | 38,480 |
| After level 2 (10µm) | 192,753 | 213,200 |
| After level 3 (5µm) | **1,424,746** | 1,610,960 |

Actual came in ~88% of the continuous-box estimate (topoSet's box selects
whole cells by center, so the true count is always a bit off from the
idealized continuous volume) — close enough that the estimates above are
trustworthy for planning.

**Mesh quality:** `checkMesh` reports **Mesh OK** — max non-orthogonality
26.1 (well within normal limits), max skewness 0.33, max aspect ratio 1.04.
`refineHexMesh` automatically inserts polyhedral cells (44,577 of them, ~3%
of the total) at the refinement-level boundaries to keep the mesh conformal
— i.e. **no manual buffer/grading layers were needed**; the "is it
straightforward" question from earlier is answered: yes, `topoSet` +
`refineHexMesh` handles this cleanly out of the box.

**Memory (mesh generation + decomposePar only, serial, cgroup `memory.peak`):**
blockMesh + 3× refine + checkMesh peaked at **~1.25 GiB**; `decomposePar`
(16-way) peaked at **~1.09 GiB**, ~14s wall time. Both are trivial next to
the old `--memory=76g` budget — but **this does not measure solver memory**
(field storage, matrix assembly, PIMPLE state across 1.4M cells is a
different and much larger cost than bare mesh topology). Still open: an
actual `laserbeamFoam` run (even a few timesteps) to get a real solver
memory reading before finalizing `--memory`/`--cpus` for the production runs.

### Second refinement box: surface overflow (2026-07-23)

Added a second static box to also resolve material overflow at the surface
(spatter/flow spreading laterally beyond the melt pool proper): y = 55µm
above the surface to 5µm below it (0.145–0.205mm, 60µm total — narrower than
an earlier 100-above/5-below draft), x = ±250µm (widened from an initial
±200µm draft), z full domain. Both boxes are unioned into one `topoSetDict`
(`new` + `add` boxToCell actions) and refined together to 5µm.

### Mesh-generation test (`testrun56_vdep_mesh_test`, 2026-07-23)

Same pipeline as testrun55, with `topoSetDict` now selecting the union of
both boxes:

| Stage | Cells |
|---|---|
| Base mesh | 16,640 |
| After level 1 | 41,574 |
| After level 2 | 239,779 |
| After level 3 (5µm, both boxes) | **1,774,977** |

`checkMesh`: **Mesh OK** (non-orthogonality max 26.1, skewness 0.33, aspect
ratio 1.04 — same healthy quality as the single-box case). Peak memory
(mesh generation only): ~1.52 GiB. At `maxCells=2,500,000` this leaves
~725k cells of headroom for `dynamicRefineFvMesh` beyond the two static
boxes — still workable, though less than testrun55's ~1.08M-cell headroom;
worth revisiting `maxCells` once a real solver-memory reading is in hand.

`report/plot_domain_schematic.py` was updated to match: clipped domain
(y 0→0.5mm), both refinement boxes drawn as dashed outlines (box 2 at the
current ±250µm), the laser track extended in black (no legend) from the
active red segment out to the domain's z edges, and the axis ticks curated
to only meaningful values (domain extent, refinement-box bounds on x; domain
extent + the three beam-event z-locations on z, which replaced the
per-event red text labels).

### Box 1 shrunk further (`testrun57_vdep_mesh_test`, 2026-07-23)

Two more trims to box 1, since box 2 already covers the near-surface/gas
region:

- **y**: no longer extends above the surface at all — now 0.2mm (surface)
  to 0.4mm (200µm below), was 0.145mm to 0.4mm.
- **z**: limited to 0.7mm–2.7mm (a 2mm window around the track), was full
  domain (0–3.2mm).

Box 2 unchanged (y 0.145–0.205mm, x ±250µm, full z).

| Stage | Cells |
|---|---|
| Base mesh | 16,640 |
| After level 1 | 32,908 |
| After level 2 | 169,002 |
| After level 3 (5µm, both boxes) | **1,203,273** |

Real, measured (same `topoSet`+`refineHexMesh`×3 pipeline). `checkMesh`:
**Mesh OK** (non-orthogonality max 26.1, skewness 0.33, aspect ratio 1.04 —
unchanged quality). Peak memory (mesh generation only): ~1.05 GiB. Down
~32% in cell count from testrun56 (1,774,977 → 1,203,273) from shrinking
and z-limiting box 1 alone.

`report/plot_domain_schematic.py` updated to match (box 1 now drawn
shorter in z, starting flush with the metal surface instead of poking
above it).

### Box 2 also z-capped at 2700µm (`testrun57_vdep_mesh_test`, 2026-07-23)

Box 2 (surface overflow) was still full-z (0–3.2mm); capped its upper bound
to match box 1's, z ≤ 2.7mm (lower bound unchanged, still unrestricted/full
from z=0). Box 2 is now y[145,205]µm, x±250µm, z ≤ 2700µm.

| Stage | Cells |
|---|---|
| Base mesh | 16,640 |
| After level 1 | 31,879 |
| After level 2 | 157,207 |
| After level 3 (5µm, both boxes) | **1,114,408** |

Real, measured. `checkMesh`: **Mesh OK** (same quality as prior meshes).
Peak memory (mesh generation only): ~0.79 GiB. Down a further ~7% from
1,203,273 → 1,114,408. Applied to all 9 case directories (testrun57–65) and
the schematic.

### `maxCells` finalized at 2,500,000 (2026-07-23)

| | Cells |
|---|---|
| Static floor (testrun56, both boxes, box 1 full-z @5µm) | 1,774,977 |
| Static floor (testrun57 pre-z-cap, box 1 z-limited only @5µm) | 1,203,273 |
| Static floor (testrun57 current, both boxes z-limited @5µm) | 1,114,408 |
| `maxCells` budget | 2,500,000 |
| **Headroom vs. current floor** | **~1,385,592 (≈55% of budget)** |

## Task 2 — Build the staged case set (2026-07-23)

Built the full 8-case staged pipeline (mirrors the tr50→tr51→tr52/53 pattern
from `report/jul13-tasks.md`, on the new pinned mesh):

| Case | Stage | Time window | Power | Mesh | Notes |
|---|---|---|---|---|---|
| `testrun58_vdep_3_Al` | 1 (tr50-equivalent) | 0–20µs | 1000 W | builds fresh: `blockMesh` → 3×(`topoSet`+`refineHexMesh`) → `setFields` → `decomposePar` → run | only case with `initial/` + full `Allrun` |
| `testrun59_vdep_3_Al` | 2 (tr51-equivalent) | 20–100µs | 650 W | continuation | seed from testrun58's last timestep |
| `testrun60_vdep_3_Al` | 3, branch | 100–400µs | 650 W | continuation | seed from testrun59 |
| `testrun61_vdep_3_Al` | 3, branch | 100–400µs | 700 W | continuation | seed from testrun59 |
| `testrun62_vdep_3_Al` | 3, branch | 100–400µs | 750 W | continuation | seed from testrun59 |
| `testrun63_vdep_3_Al` | 3, branch | 100–400µs | 800 W | continuation | seed from testrun59 |
| `testrun64_vdep_3_Al` | 3, branch | 100–400µs | 850 W | continuation | seed from testrun59 |
| `testrun65_vdep_3_Al` | 3, branch | 100–400µs | 900 W | continuation | seed from testrun59 |

All 8 share the pinned geometry/mesh (clipped 500µm domain, both static
refinement boxes per `testrun57`'s bounds, `dynamicMeshDict` with
`maxCells=2,500,000`, `field=alpha.metal` — the unrefine risk noted below is
still unresolved, carried into every case). `testrun58`'s `Allrun` runs the
static refinement loop directly (not via `runApplication`, which would skip
passes 2–3 since the log file from pass 1 already exists). The 7
continuation cases only differ in `constant/timeVsLaserPower` and
`system/controlDict`'s `endTime`; per TESTRUNS.md's documented multi-stage
pattern, none of them have been run yet — each needs its predecessor's
final timestep manually copied in first (testrun58 must actually be run
before testrun59 can be seeded, etc.).

**Still open / not yet done:**
- Get a real solver-memory reading (a short `laserbeamFoam` run) before
  fixing `--memory`/`--cpus` for the actual runs.
- The combined static-region + `alpha.metal` trigger field for
  `dynamicRefineFvMesh` (box 2's gas-side unrefine risk, below) — not
  implemented; all 8 cases currently use plain `field alpha.metal`.

**Open technical risk (unresolved, carried into all 8 cases):**
`dynamicRefineFvMesh`'s unrefine rule only fires when
`alpha.metal < unrefineLevel` (0.001) — i.e. it only coarsens cells that
read as pure gas. Bulk-metal cells (α=1, most of box 1's depth) are never
marked for unrefinement so they'll stay fine indefinitely, which is
fine/intended. But gas-side cells inside box 2 (the surface-overflow box,
mostly α≈0 away from the melt pool) *do* satisfy `α < 0.001` and would be
marked for unrefinement over time, undoing the static pre-refinement there
as the run progresses — unless the trigger field is changed to a combined
static-region + `alpha.metal` indicator so the solver knows those cells are
protected regardless of local α.

### Box 2 z-capped to match box 1 (2026-07-23)

Box 2 (surface overflow) was still full-z; capped its upper bound to match
box 1's, z ≤ 2.7mm (lower bound unchanged, still unrestricted from z=0).
Applied to all 9 case directories (testrun57–65) and the schematic. Real
mesh (`testrun57`):

| Stage | Cells |
|---|---|
| Base mesh | 16,640 |
| After level 1 | 31,879 |
| After level 2 | 157,207 |
| After level 3 (5µm, both boxes) | **1,114,408** |

`checkMesh`: **Mesh OK**. Peak memory (mesh generation only): ~0.79 GiB.
Headroom vs. `maxCells=2,500,000`: **~1,385,592 (≈55%)** — the best margin
of any iteration.

## Task 3 — Solver fixes, recording fields, and output format (2026-07-23)

### Gas-branch kappa now also clamped at Tvap

`applications/solvers/laserbeamFoam/updateKappaCp.H` previously only
clamped the *metal liquid* branch's kappa-polynomial argument at `Tvap`
(commit `223825ed`, 2026-07-19) — the gas branch still evaluated
`polykappa_g` with raw, unclamped `T`. Extended the same clamp to the gas
branch so both use `TKappaClamped = min(T, Tvap)`:

```cpp
const scalar TKappaClamped = min(TI[celli], TvapValue);
...
+ epsilon1I[celli]*polykappa_m_liquidPtr().value(TKappaClamped)   // metal liquid branch
...
(1.0 - alpha1I[celli]) * polykappa_g.value(TKappaClamped)         // gas branch (was raw T)
```

No numerical effect on the current cases today — `transportProperties`'
gas `poly_kappa` is a flat constant `(0.04 0 0 0 0 0 0 0)`, no T-dependence
— but protects against the same runaway/instability risk if a
T-dependent gas kappa is ever used.

### Field-write suppression (reduces recording size)

No dict-level mechanism exists for suppressing individual `AUTO_WRITE`
fields in this solver (`laserbeamFoam.C` just calls plain `runTime.write()`,
which writes every registered `AUTO_WRITE` field — confirmed by checking
for any `writeObjects`/allowlist precedent across the repo's tutorials;
none exists). Flipped 5 fields from `AUTO_WRITE` to `NO_WRITE` directly in
`createFields.H` instead — all confirmed `NO_READ` (never used to resume a
run, so this doesn't affect the staged-continuation seeding):

| Field | Why safe to drop |
|---|---|
| `TSolidus`, `TLiquidus` | Static per-cell values (alpha1-blended constants) |
| `cp`, `kappa` | Recomputed from `T`/`alpha.metal`/`epsilon1` every step |
| `Qv` | Recomputed every step |

Remaining recorded fields per snapshot: `alpha.metal`, `U`, `T`, `p_rgh`,
`p`, `epsilon1`, `alpha_smoothed`, `alpha_filtered`, `meltHistory`,
`condition` (10, down from 15).

### Output format: binary + compression

All 8 case `controlDict`s switched from `writeFormat ascii` /
`writeCompression off` to **`binary`** / **`on`**. ParaView's OpenFOAM
reader handles both natively — no downside for reading, and it's the
standard combination for minimizing disk footprint (binary avoids ASCII's
per-number text overhead; gzip compresses further on top).

### Solver rebuilt

`docker build --build-arg CACHE_BUST=... -t lbf3 .` run twice (once after
the gas-kappa clamp, once after the field-suppression edit) — both clean,
**"There were no build errors"** (`laserbeamFoam`, `compressibleLaserbeamFoam`,
`multicomponentLaserbeamFoam`, `laserMeltFoam`, `setSolidFraction` all
compiled). Current image must be used for any run of the 8 cases — an
older `lbf3` image would be missing both fixes.

**Recording schedule** (all 8 cases, `writeInterval=4µs`, `purgeWrite=0`):

| Stage | Snapshots in window | Cumulative |
|---|---|---|
| Seed 1 (testrun58, 0–20µs) | 6 (0, 4, …, **20**) | 6 |
| Seed 2 (testrun59, 20–100µs) | 20 new (24, …, **100**) | 26 |
| Each branch (testrun60–65, 100–400µs) | 75 new each (104, …, **400**) | **101/branch** |

The bolded stage-boundary times are guaranteed on disk (OpenFOAM always
writes at `endTime`), which is what makes them valid seeds for the next
stage regardless of `writeInterval` alignment.

**Still open:**
- Real solver-memory reading (no case has been run yet).
- Combined static-region + `alpha.metal` trigger field for
  `dynamicRefineFvMesh` (box 2 gas-side unrefine risk, above) — still
  unresolved.
