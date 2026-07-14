# Jul 13 Tasks — metal kappa(T) fix + testrun50/51/52/53 staged runs

## Context

`testrun50`'s metal-phase thermal conductivity (`poly_kappa`) was a single global
degree-7 polynomial forced to cover the whole temperature range: flat ~167 W/(m·K)
below the solidus, then a physical drop as the metal melts, then a gentle decline
in the liquid up to ~2743 K. We proved via minimax (Chebyshev) optimization that no
single degree-7 polynomial can represent that shape with better than ~25.7 W/(m·K)
worst-case error — a hard mathematical floor, not a fitting mistake. The fix is to
stop asking one polynomial to do a job it structurally can't, and instead blend a
solid-branch and a liquid-branch polynomial using the solver's own liquid-fraction
field (`epsilon1`), which is already computed for the mushy-zone Darcy damping.

We're also splitting testrun50's original single 0–100 µs run into a staged
pipeline so each stage's melt-pool physics can be checked independently before
continuing, and forking the final stage across two power levels for comparison —
run by two people in parallel:

- **testrun50** (stage 1): 0 → 20 µs, constant 1000 W — run separately by both Mehrdad and Zixun
- **testrun51** (stage 2): 20 → 100 µs, constant 650 W, continues from each person's own testrun50 — run separately by both
- **testrun52** (stage 3, Zixun's branch): 100 → 400 µs, constant **800 W**, continues from Zixun's testrun51
- **testrun53** (stage 3, Mehrdad's branch): 100 → 400 µs, constant **850 W**, continues from Mehrdad's testrun51

## Tasks

- [x] **Fix metal kappa(T) to be physically accurate (Option 4: piecewise blend)**
  - Added optional `poly_kappa_liquid` entry to `createFields.H`, guarded by
    `transportPropertiesMetal.found(...)` (mirrors the existing `poly_sigma`
    pattern) — fully backward-compatible, all other tutorial cases lacking
    this key take the old single-polynomial code path unchanged.
  - `updateKappaCp.H` now blends solid/liquid kappa by the local liquid
    fraction `epsilon1` when the liquid entry is present.
  - Reordered `TSolidus`/`TLiquidus`/`epsilon1` declarations in `createFields.H`
    ahead of `cp`/`kappa` — a compile error (`epsilon1` not yet declared at the
    first `updateKappaCp.H` include) caught this ordering issue immediately;
    fixed and rebuilt clean.
  - `testrun50/constant/transportProperties`: `poly_kappa` is now flat solid
    `(167 0 0 0 0 0 0 0)`; new `poly_kappa_liquid` fitted to 95 W/(m·K) at
    Tliquidus declining to ~78 W/(m·K) near Tvap (minimax fit error ~8e-5
    W/(m·K) — essentially exact).
  - Verified: clean Docker rebuild, and confirmed untouched cases (testrun44,
    testrun1) still lack `poly_kappa_liquid` and take the old code path.
  - Verified in a live run: testrun50 ran for >1 hour with converging
    residuals, bounded Courant numbers, and a monotonically converging
    `epsilon1` (liquid fraction) correction each PIMPLE iteration — the exact
    field the new kappa blend depends on. No instability observed.

- [x] **testrun50: restrict to stage 1 (0–20 µs, constant 1000 W)**
  - `system/controlDict`: `endTime` 1e-4 → 2e-5.
  - `constant/timeVsLaserPower`: two-stage (1000 W → 500 W at 20 µs) →
    constant 1000 W for the whole stage.
  - Laser position table (same 400 µs / 2.4 mm track, 6 m/s) left unchanged —
    shared across all stages so the beam continues seamlessly.
  - Currently running (started by Mehrdad).

- [x] **testrun51: stage 2 (20–100 µs, constant 650 W), continuing from testrun50**
  - New case directory, copied from testrun50's structure (same mesh,
    material, AMR, decomposition).
  - `constant/timeVsLaserPower` set to constant 650 W.
  - `system/controlDict`: `endTime` = 1e-4; `startFrom latestTime` (already
    the case default) picks up whatever time directory is copied in.
  - `Allrun` stripped down to a continuation-only script: no `initial → 0`
    copy, no `blockMesh`, no `setFields`, no `decomposePar` — just runs the
    solver directly against data copied forward from testrun50's final
    timestep.
  - Requires testrun50 to actually be run first so there's a 20 µs checkpoint
    to copy forward. Each person runs their own copy of testrun50 → testrun51
    (not shared — this becomes the seed for their own branch below).

- [x] **testrun52: stage 3, Zixun's branch (100–400 µs, constant 800 W), continuing from testrun51**
  - Same setup pattern as testrun51 (continuation-only `Allrun`).
  - `constant/timeVsLaserPower` set to constant 800 W (changed from an
    earlier draft at 650 W).
  - `system/controlDict`: `endTime` = 4e-4.
  - Requires testrun51 to be run first so there's a 100 µs checkpoint to copy
    forward. Run by Zixun.

- [x] **testrun53: stage 3, Mehrdad's branch (100–400 µs, constant 850 W), continuing from testrun51**
  - Same setup as testrun52, copied from it, with power changed to 850 W.
  - Continuation-only `Allrun` (inherited from testrun52).
  - `system/controlDict`: `endTime` = 4e-4.
  - Run by Mehrdad, continuing from his own testrun51 result.
