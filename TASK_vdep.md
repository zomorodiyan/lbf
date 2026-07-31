# Task: Replicate VDEP Experimental Results (AlSi10Mg, laserbeamFoam)

## Ultimate Goal

Reproduce experimentally measured vapor depression (keyhole) depth and morphology
for AlSi10Mg single-track laser scans across a range of powers, using the
`laserbeamFoam` VOF solver. Validation target: keyhole depth vs. laser power curve
from X-ray synchrotron or ex-situ characterization data.

---

## Current Pipeline

The active case set is **testrun58–testrun67** — a two-stage seed→fork lineage (a
shared 100µs warmup/stabilization state at each of two core counts, then six power
forks branching off it). Full case table, per-person assignments, and the exact
seed→fork hand-off procedure live in [vdep_power_sweep.md](vdep_power_sweep.md) —
read that before running, seeding, or reconstructing any of these cases. This file
covers the physical setup, overall goal, and lessons that apply across the whole
effort.

## Physical Setup

- **Material:** AlSi10Mg (T_sol=840 K, T_liq=867 K, ρ=2670 kg/m³)
- **Laser:** 1064 nm, 35 µm radius, Drude/Fresnel absorptivity (~11% at normal
  incidence for liquid Al)
- **Scan speed:** 6000 mm/s (6 m/s)
- **Powers of interest:** 650–900 W (six forks, see case table in
  vdep_power_sweep.md)
- **Mesh:** 40 µm base cells, 3-level AMR at the metal–gas interface → ~5 µm
- **Decomposition:** hierarchical — 16 cores (n=(2 1 8)) or 32 cores (n=(2 2 8)),
  depending on the fork; must stay fixed across an entire seed→fork lineage

---

## Progress

- [x] Seed lineage built and validated: testrun58/59 (16-core) and testrun60/61
  (32-core), each running 0→100µs (1000W warmup, then 650W stabilization)
- [x] Seed→fork mesh/field hand-off procedure tested end-to-end
  (testrun58→testrun59) — see vdep_power_sweep.md's "Hand-off validated" note
- [ ] Reconstruct and post-process all six forks once complete
- [ ] Extract keyhole depth vs. power and compare to experiment

### Fork status (as of 2026-07-31)

| Case | Power | Owner | Progress | Notes |
|---|---|---|---|---|
| testrun62 | 650 W | Mehrdad | 80% (t≈340µs) | running |
| testrun63 | 700 W | Mehrdad | 0% | forked, not yet started |
| testrun64 | 750 W | Mehrdad | 100% (t=400µs) | complete |
| testrun66 | 850 W | Mehrdad | 15% (t≈145µs) | running |
| testrun65 | 800 W | Zixun | 73% (t≈320µs) | running, latest frame reconstructed |
| testrun67 | 900 W | Zixun | 60% (t≈280µs) | running, latest frame reconstructed |

Zixun is relocating and packing up his PC, so his animation rendering for
testrun65/67 may be delayed until he's set back up — reconstructed data is
ready whenever he (or anyone else) gets to it.

Check with whoever owns a given fork for anything more recent than this table
— it's a snapshot from the messages above, not a live status.

---

## Post-processing: turning reconstructed results into animations

`results/_render_stacked_video.sh` batch-renders every reconstructed timestep of a
case into three views (top-down, lateral, synthetic X-ray), stacks them into one
image per timestep, and assembles an mp4. Full details, including the exact script
options and troubleshooting, are in TESTRUNS.md's "Post-processing" and
"Permissions across both Docker images" sections — short version:

1. Make sure the case is reconstructed (`reconstructParMesh` + `reconstructPar` —
   see TESTRUNS.md).
2. One-time setup, per machine:
   - Rebuild `lbf3` if it predates the `ffmpeg` addition to the Dockerfile:
     `docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .`
   - Pull the second image the script needs (public, not something we build):
     `docker pull kitware/paraview:pv-v5.8.0-osmesa-py3`
3. Run it: `bash results/_render_stacked_video.sh <testrun-number-or-path>`, e.g.
   `bash results/_render_stacked_video.sh 65`.

If you hit a `Permission denied` error partway through, it's almost always a
Docker user-mismatch issue, not a real problem with the data — see TESTRUNS.md's
"Permissions across both Docker images" section for the fix.

---

## Key Lessons Learned

- Running multiple cases simultaneously on one machine causes severe slowdown
  (~12x observed with 4 concurrent runs) due to memory bandwidth contention —
  run sequentially, or at most 2 at a time.
- `processor*/` directories are root-owned (written by Docker as root) — always
  delete them via Docker, never `sudo rm` (see TESTRUNS.md).
- Decomposition (core count) must stay identical across an entire seed→fork
  lineage — continuation `Allrun` scripts read their own `decomposeParDict`
  rather than inspecting the `processor*/` dirs they were seeded with, so a
  core-count change can only happen at a fresh seed stage, not mid-lineage.
