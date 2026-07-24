# VDEP power-sweep pipeline (testrun58–testrun67)

How to run the active VDEP power-sweep case set, including Zixun's portion. Geometry/mesh
derivation: `report/jul23-tasks.md`. Generic Docker mechanics: [TESTRUNS.md](TESTRUNS.md).

## Case table

| Case | Role | Time window | Power | Cores | Seeded from |
|---|---|---|---|---|---|
| `testrun58_vdep_3_Al` | seed-0_16c | 0–20µs | 1000 W | 16 | fresh build |
| `testrun59_vdep_3_Al` | seed-1_16c | 20–100µs | 650 W | 16 | testrun58 |
| `testrun60_vdep_3_Al` | seed-0_32c | 0–20µs | 1000 W | 32 | fresh build |
| `testrun61_vdep_3_Al` | seed-1_32c | 20–100µs | 650 W | 32 | testrun60 |
| `testrun62_vdep_3_Al` | fork | 100–400µs | **650 W** | 16 | testrun59 |
| `testrun63_vdep_3_Al` | fork | 100–400µs | **700 W** | 16 | testrun59 |
| `testrun64_vdep_3_Al` | fork | 100–400µs | **750 W** | 32 | testrun61 |
| `testrun65_vdep_3_Al` | fork | 100–400µs | **800 W** | 32 | testrun61 |
| `testrun66_vdep_3_Al` | fork | 100–400µs | **850 W** | 32 | testrun61 |
| `testrun67_vdep_3_Al` | fork | 100–400µs | **900 W** | 32 | testrun61 |

The six bolded forks (testrun62–67) are the actual production runs — the seed cases (58–61)
only exist to produce the shared 100µs starting point each fork branches from, at the two core
counts the forks need.

## Assignment

- **Mehrdad** — 650W (testrun62) and 700W (testrun63) on the 20-core workstation (16 cores used),
  seeded from testrun59 (16-core lineage). 750W (testrun64) and 850W (testrun66) on the
  blackwell workstation (32 cores), seeded from testrun61 (32-core lineage).
- **Zixun** — 800W (testrun65) and 900W (testrun67), seeded from testrun61 (32-core lineage).
  Needs a 32-core-capable machine to match that lineage's decomposition.

## Zixun's steps

1. Build and run **testrun60** (seed-0_32c) — a fresh build, same as any seed case:
   ```bash
   cd ~/lbf3
   docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
     --ipc=host --cpus=32 --memory=76g \
     -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun60_vdep_3_Al && bash ./Allrun && echo RUN_COMPLETE"
   ```
   Tail with `tail -f tutorials/laserbeamFoam/vdep/testrun60_vdep_3_Al/log.laserbeamFoam`.
   Runs to `endTime=20µs`.

2. **Copy testrun60's final state into testrun61.** Continuation cases don't have `constant/polyMesh`
   or any decomposed `processor*/` of their own yet — both come from the seed. Procedure:
   ```bash
   # a. Reconstruct testrun60's latest timestep
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun60_vdep_3_Al && \
      reconstructParMesh -noFunctionObjects -latestTime 2>&1 | tee log.reconstructParMesh && \
      reconstructPar    -noFunctionObjects -latestTime 2>&1 | tee log.reconstructPar"

   # b. Note the reconstructed time directory name (e.g. 2e-05), then copy the base mesh
   #    and that one timestep into testrun61
   cp -r tutorials/laserbeamFoam/vdep/testrun60_vdep_3_Al/constant/polyMesh \
         tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al/constant/
   cp -r tutorials/laserbeamFoam/vdep/testrun60_vdep_3_Al/<TIME> \
         tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al/<TIME>

   # c. Decompose that timestep in testrun61 using ITS OWN decomposeParDict (32-core)
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al && decomposePar -time <TIME>"

   # d. Sanity-check the copied/decomposed mesh before committing hours of compute to it
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al && checkMesh -time <TIME>"
   ```
   Expect `Mesh OK` with the same domain bounding box, non-orthogonality (~26), and skewness
   (~0.33) as the source case. Don't proceed to step 3 if this doesn't come back clean.

3. Run **testrun61** (seed-1_32c) — now that `processor*/` exists, its `Allrun` will skip
   straight to `runParallel`:
   ```bash
   docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
     --ipc=host --cpus=32 --memory=76g \
     -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al && bash ./Allrun && echo RUN_COMPLETE"
   ```
   Runs 20µs→100µs at 650W.

4. Repeat the same reconstruct → copy `constant/polyMesh` + latest timestep → `decomposePar`
   hand-off (step 2) from **testrun61 into testrun65**, then run testrun65's `Allrun`
   (100µs→400µs at 800W).

5. Repeat the identical hand-off from **testrun61 into testrun67**, then run testrun67's `Allrun`
   (100µs→400µs at 900W).

   testrun65 and testrun67 are independent forks of the same testrun61 result — the hand-off in
   step 4 and step 5 both start from testrun61, not from each other.

## Hand-off validated (2026-07-24)

The reconstruct→copy `constant/polyMesh`+timestep→`decomposePar` procedure above was tested
end-to-end on testrun58→testrun59: `reconstructParMesh`/`reconstructPar -latestTime` on testrun58
(16 processors, t=2e-05), copy into testrun59, `decomposePar -time 2e-05` with testrun59's own
16-core `decomposeParDict`, then `checkMesh` (Mesh OK — same domain bounding box, non-orthogonality
26.1, skewness 0.33 as the source) before launching. testrun59 then ran correctly: laser power
switched to 650W as expected and `alpha.metal` phase fraction continued smoothly from testrun58's
end state (≈0.603) rather than resetting — confirming the mesh/field hand-off actually carries the
run forward rather than silently restarting from a wrong state. Use the same steps with confidence
for testrun60→61 and testrun61→65/67.
