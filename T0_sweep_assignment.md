# T0 sweep run assignment (testrun73–testrun84)

Who runs what, and in what order, for the four T0 lineages prepared in
[T0_sweep.md](T0_sweep.md) — read that doc first for what these cases are, why the fork stage
carries testrun69's surface-tension remedy, and the exact reconstruct→copy→decompose hand-off
procedure (same one used in [vdep_power_sweep.md](vdep_power_sweep.md)). This doc only covers
who runs which lineage and in what sequence.

## Split

- **Mehrdad** — cold pair, coldest first:
  1. T0=100K: `testrun73` → `testrun74` → `testrun75`
  2. T0=200K: `testrun76` → `testrun77` → `testrun78`
- **Zixun** — hot pair, hottest first:
  1. T0=500K: `testrun82` → `testrun83` → `testrun84`
  2. T0=400K: `testrun79` → `testrun80` → `testrun81`

## Ordering rules

- Within one T0 lineage, the three stages are strictly sequential: seed0 must finish and be
  reconstructed before its seed1 hand-off can happen, and seed1 must finish and be reconstructed
  before the fork hand-off. You cannot skip ahead or run them out of order.
- Across your own two T0 lineages, finish your first one (coldest for Mehrdad, hottest for Zixun)
  before starting the second — per TASK_vdep.md's "Key Lessons," running multiple cases at once on
  one machine causes severe slowdown (~12x observed with 4 concurrent), so don't run both of your
  lineages' stages in parallel on the same box.
- All twelve cases use the same 32-core decomposition as testrun61's lineage — same
  `--cpus=32 --memory=...` Docker flags as any other 32-core VDEP case (see TESTRUNS.md).

## Per-stage procedure (same steps, three times per lineage)

**Stage 1 (seed0 — e.g. testrun73)** — fresh build, runs like testrun60:
```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al && bash ./Allrun && echo RUN_COMPLETE"
```
Runs to `endTime=20µs`.

**Hand-off (seed0 → seed1, or seed1 → fork)** — identical reconstruct→copy→decompose→checkMesh
procedure documented in vdep_power_sweep.md's "Hand-off validated" section. Reconstruct the source
stage's latest timestep, copy `constant/polyMesh` + that timestep into the next stage, decompose
with the next stage's own `decomposeParDict`, `checkMesh` before committing compute. **Seed from
your own T0's predecessor stage, not from testrun61/testrun69** — e.g. testrun74 seeds from
testrun73, testrun84 seeds from testrun83.

**Stage 2 (seed1) / Stage 3 (fork)** — once seeded, `Allrun` skips straight to `runParallel`:
```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al && bash ./Allrun && echo RUN_COMPLETE"
```
seed1 runs to `endTime=100µs` (650W); the fork runs to `endTime=400µs` (750W + sigma remedy).

## After each lineage finishes

Reconstruct (`reconstructParMesh` then `reconstructPar`, see TESTRUNS.md) and post-process with
`results/_render_stacked_video.sh <NN>`, same as the power-sweep forks.
