# VDEP recoil-pressure remedy sims (testrun68-70)

How to run the three remedy simulations from [task_Aug4.md](task_Aug4.md), split between Mehrdad
and Zixun. Physics/parameter rationale, literature grounding, and supporting plots live in
task_Aug4.md — this file only covers the run mechanics. Generic Docker mechanics:
[TESTRUNS.md](TESTRUNS.md). Seed→fork hand-off precedent (same procedure reused below):
[vdep_power_sweep.md](vdep_power_sweep.md).

## Case table

| Case | Remedy | Change | Owner |
|---|---|---|---|
| `testrun68_vdep_3_Al` | Recoil pressure cap | `Tcap = 3800 K` (new case-configurable parameter) | Mehrdad |
| `testrun69_vdep_3_Al` | Surface tension | `sigma: 0.87 → 0.95 N/m` | Mehrdad |
| `testrun70_vdep_3_Al` | Beam quality | `Radius_Flavour: 1.336 → 0.7` | Zixun |

All three fork from **testrun61** (100µs, 32-core seed — same lineage as the production forks
testrun64-67) and are otherwise **identical to testrun64's configuration** (750W,
`endTime`=400µs, same decomposition, same everything else) — only the one remedy parameter
differs per case. This isolates each remedy's effect against testrun64 as a direct, real
baseline, per task_Aug4.md's "one remedy per sim" design.

---

## ⚠ testrun68 needs its own Docker image — Mehrdad-only, nothing for Zixun to do

The recoil-pressure cap is a **solver-code change**, not a case-dictionary edit like the other
two remedies. This is entirely **Mehrdad's side** — Zixun never touches this section, never
rebuilds anything, and keeps using the plain `lbf3` image exactly as before for testrun70.

**Build it as a separate, distinctly-named image — do NOT rebuild the shared `lbf3` image in
place**:
```bash
docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3-tcap .
```
(after making the code change below on a branch/copy that has the `Tcap` clamp added). Keeping it
as `lbf3-tcap` rather than overwriting `lbf3` means nobody else's in-progress work on the shared
image is disturbed, and there's no coordination needed with Zixun before doing this.

**Reminder for later**: every `docker run` command for testrun68 (build, resume, reconstruct —
anything touching that case) must say `lbf3-tcap`, not `lbf3`. Everything else in this file
(testrun69, testrun70, and the seed hand-off steps below for all three cases) still uses the
regular `lbf3` image — only testrun68's *own* run/resume commands swap the image name.

Code change itself (`Tcap` as a case-configurable dictionary parameter, looked up the same way
`Tvap` already is, defaulting to "no cap" when absent — so even reused on `lbf3-tcap`, a case
without a `Tcap` entry behaves exactly like the old unbounded formula): add a `Tcap` lookup next
to `Tvap`'s in `createFields.H`, then clamp `TsafeU` in `UEqn.H` at both the existing 300K floor
and the new `Tcap` ceiling before it feeds the recoil-pressure exponent. Exact lookup syntax TBD
when actually coding it.

---

## Per-case setup

Repeat for each of testrun68/69/70:

1. **Copy testrun64 as the template** (identical config — `controlDict`, `decomposeParDict`,
   `fvSchemes`, `fvSolution`, `transportProperties`, `LaserProperties`, `topoSetDict`,
   `timeVsLaserPower` @750W):
   ```bash
   cp -r tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al
   ```

2. **Strip testrun64's own run state** from the copy — the new case must reseed from testrun61's
   100µs state, not continue from testrun64's already-evolved 400µs end state:
   ```bash
   cd tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al
   ls    # check first -- confirm what's actually there before deleting anything
   rm -rf processor* [0-9.e-]* log.* VTKs *.foam
   # constant/ and system/ must survive -- that's the template config being preserved
   ```

3. **Seed hand-off from testrun61** — identical procedure already validated for
   testrun58→59 and reused for testrun61→64/65/66/67 (see
   [vdep_power_sweep.md](vdep_power_sweep.md)'s "Hand-off validated" section):
   ```bash
   # a. Reconstruct testrun61's latest timestep (skip if already done)
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al && \
      reconstructParMesh -noFunctionObjects -latestTime && reconstructPar -noFunctionObjects -latestTime"

   # b. Copy the mesh + that one timestep into testrunNN
   cp -r tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al/constant/polyMesh \
         tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al/constant/
   cp -r tutorials/laserbeamFoam/vdep/testrun61_vdep_3_Al/<TIME> \
         tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al/<TIME>

   # c. Decompose using testrunNN's own (copied-from-tr64) decomposeParDict
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al && decomposePar -time <TIME>"

   # d. Sanity check before committing compute
   docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
     "cd /workspace/tutorials/laserbeamFoam/vdep/testrunNN_vdep_3_Al && checkMesh -time <TIME>"
   ```
   Expect `Mesh OK`, matching testrun61's own bounding box/non-orthogonality/skewness — don't
   proceed if this doesn't come back clean.

4. **Add the mass-conservation monitor** to `system/controlDict`'s `functions` block (not present
   in testrun64 — new, per task_Aug4.md's Tasks 2-4 setup):
   ```
   functions
   {
       parProfiling { type parProfiling; libs ( utilityFunctionObjects ); }
       metalVolume
       {
           type          volFieldValue;
           libs          ( fieldFunctionObjects );
           fields        ( alpha.metal );
           operation     volIntegrate;
           regionType    all;
           writeFields   false;
           writeControl  writeTime;
       }
   }
   ```
   `writeFields` is mandatory for this function object in OpenFOAM v2506 — omitting it is a
   `FOAM FATAL IO ERROR` at the very start of the run (confirmed the hard way on testrun69).

5. **Apply the one remedy-specific change** — see below.

6. **Run** — same command as any other VDEP case, see [TESTRUNS.md](TESTRUNS.md)'s "Run" section.

---

## Per-case remedy change

### testrun68 — recoil pressure cap
Requires the `lbf3-tcap` image above (Mehrdad-only — see the ⚠ section). Add to
`constant/transportProperties`:
```
Tcap 3800;
```
Run/resume/reconstruct commands for this case use `lbf3-tcap` in place of `lbf3` — everything
else about the command is unchanged.

### testrun69 — surface tension
No solver change needed — runs on the stock `lbf3` image. Edit `constant/transportProperties`:
```
sigma 0.95;  // was 0.87
```

### testrun70 — beam quality
No solver change needed. Edit `constant/LaserProperties`:
```
Radius_Flavour 0.7;  // was 1.336
```

---

## Assignment

- **Mehrdad** — testrun68 (recoil pressure cap — including the `lbf3-tcap` image build, entirely
  on this side) and testrun69 (surface tension).
- **Zixun** — testrun70 (beam quality / Radius_Flavour) only. Uses the plain `lbf3` image, same
  as every other case — nothing solver-related to build or coordinate on. Needs the
  32-core-capable machine used before for testrun65/67 (same decomposition as testrun61's
  lineage).
