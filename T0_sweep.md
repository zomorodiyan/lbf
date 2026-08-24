# VDEP initial-temperature sweep (testrun73–testrun84)

How to run the T0-sweep case set. Generic Docker mechanics: [TESTRUNS.md](TESTRUNS.md). Seed→fork
hand-off procedure (reused verbatim below): [vdep_power_sweep.md](vdep_power_sweep.md). Remedy
baseline this sweep forks from: [vdep_remedy_sims.md](vdep_remedy_sims.md).

## Purpose

Everything in the testrun58–72 pipeline starts from ambient T0 = 300 K. This sweep repeats the
same seed→fork lineage at four other initial/reference temperatures — 100 K, 200 K, 400 K, 500 K —
to see how the vapor-depression dynamics depend on starting temperature. Each T0 gets its own
independent 3-stage lineage, structured exactly like testrun60→61→(750 W fork):

| Stage | Role | Time window | Power | Template | Notes |
|---|---|---|---|---|---|
| seed0 | fresh build | 0–20µs | 1000 W | testrun60 | own `initial/T` + `blockMesh`/`setFields` |
| seed1 | continuation | 20–100µs | 650 W | testrun61 | seeded from this T0's own seed0 |
| fork | continuation | 100–400µs | 750 W | testrun69 | seeded from this T0's own seed1 |

The fork stage is a copy of **testrun69** (not testrun64) — same 750 W power, same
`metalVolume` mass-conservation monitor. testrun69's surface-tension remedy (`sigma: 0.87→0.95`)
is applied to **every stage of every lineage** here, including seed0 and seed1 (which otherwise
came from the testrun60/61 templates at the unmodified `sigma=0.87`) — so sigma=0.95 is constant
across all three stages and all four T0 values, and only T0 itself varies step-to-step within a
lineage.

## Case table

| Case | T0 | Role | Time window | Power | Sigma | Seeded from |
|---|---|---|---|---|---|---|
| `testrun73_vdep_3_Al` | 100 K | seed0 | 0–20µs | 1000 W | 0.95 | fresh build |
| `testrun74_vdep_3_Al` | 100 K | seed1 | 20–100µs | 650 W | 0.95 | testrun73 |
| `testrun75_vdep_3_Al` | 100 K | fork | 100–400µs | 750 W | 0.95 | testrun74 |
| `testrun76_vdep_3_Al` | 200 K | seed0 | 0–20µs | 1000 W | 0.95 | fresh build |
| `testrun77_vdep_3_Al` | 200 K | seed1 | 20–100µs | 650 W | 0.95 | testrun76 |
| `testrun78_vdep_3_Al` | 200 K | fork | 100–400µs | 750 W | 0.95 | testrun77 |
| `testrun79_vdep_3_Al` | 400 K | seed0 | 0–20µs | 1000 W | 0.95 | fresh build |
| `testrun80_vdep_3_Al` | 400 K | seed1 | 20–100µs | 650 W | 0.95 | testrun79 |
| `testrun81_vdep_3_Al` | 400 K | fork | 100–400µs | 750 W | 0.95 | testrun80 |
| `testrun82_vdep_3_Al` | 500 K | seed0 | 0–20µs | 1000 W | 0.95 | fresh build |
| `testrun83_vdep_3_Al` | 500 K | seed1 | 20–100µs | 650 W | 0.95 | testrun82 |
| `testrun84_vdep_3_Al` | 500 K | fork | 100–400µs | 750 W | 0.95 | testrun83 |

All twelve use the same 32-core hierarchical decomposition as the testrun60/61 lineage (must stay
fixed across each T0's own seed→fork chain, same rule as the power sweep).

## What actually changed per case

Relative to the testrun60/61/69 templates, three things were edited:

1. **seed0 only** — `initial/T`'s `internalField   uniform 300.0;` → this T0's value. This is the
   actual initial condition; seed1 and fork stages have no `initial/` directory of their own (per
   vdep_power_sweep.md, continuation cases inherit their field state from the hand-off, not from
   a local `initial/`), so this only needs setting once per lineage, at seed0.
2. **All three stages** — `constant/transportProperties`'s `TRef` (the Boussinesq buoyancy
   reference temperature) → this T0's value. Changed in seed0, seed1, *and* the fork, since each
   stage carries its own `transportProperties`.
3. **seed0 and seed1 only** — `constant/transportProperties`'s `sigma` bumped from the
   testrun60/61 templates' `0.87` to `0.95`, to match the fork stage's testrun69-derived value.
   The fork stage already had `sigma=0.95` from its testrun69 template, so no edit was needed
   there. Net effect: sigma is `0.95` everywhere in every lineage.

Nothing else (LaserProperties, timeVsLaserPosition, timeVsLaserPower, mesh dicts, decomposeParDict)
depends on T0 — confirmed by diffing testrun73/75/84 against their testrun60/61/69 templates.

## Current status (as of 2026-08-23)

**Only the case directories/dictionaries have been prepared — nothing has been run.** Per
[CLAUDE.md](CLAUDE.md), never run a simulation without asking the user first.

- seed0 stages (testrun73/76/79/82) are complete, ready-to-run fresh builds — same
  `blockMesh` → `topoSet`×3/`refineHexMesh` → `setFields` → `decomposePar` → `laserbeamFoam`
  `Allrun` as testrun60.
- seed1 and fork stages (testrun74/75/77/78/80/81/83/84) have only their `constant/`+`system/`
  dictionaries and `Allrun`/`Allclean` copied in — they have **no `constant/polyMesh`, no
  decomposed `processor*/`, no seed timestep yet**. Each needs the same reconstruct → copy
  `constant/polyMesh` + latest timestep → `decomposePar` → `checkMesh` hand-off documented in
  vdep_power_sweep.md's "Hand-off validated" section and reused in vdep_remedy_sims.md, run from
  its own T0 lineage predecessor (e.g. testrun74 seeds from testrun73's reconstructed 20µs state,
  *not* from testrun61).

## Next steps for whoever picks this up

1. Run the four seed0 stages (testrun73/76/79/82) — can run sequentially or up to 2 at a time per
   TASK_vdep.md's concurrency guidance (severe slowdown beyond that).
2. Hand off each seed0 → its own seed1 (testrun73→74, 76→77, 79→80, 82→83), run seed1.
3. Hand off each seed1 → its own fork (testrun74→75, 77→78, 80→81, 83→84), run the fork to
   `endTime=400µs`.
4. Reconstruct and post-process with `results/_render_stacked_video.sh <NN>`, same as the power
   sweep.

Ask the user (Mehrdad) before launching any of the above — this doc only covers what was prepared,
not authorization to run.
