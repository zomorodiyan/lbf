# Claude Code instructions for this repo

## When to read TESTRUNS.md / TESTRUNS_MULTICOMPONENT.md

**Always read [TESTRUNS.md](TESTRUNS.md) (single-material `laserbeamFoam`/`laserMeltFoam`) or
[TESTRUNS_MULTICOMPONENT.md](TESTRUNS_MULTICOMPONENT.md) (`multicomponentLaserbeamFoam`,
i.e. `tutorials/multiComponentlaserbeamFoam/` and `tutorials/compressiblelaserbeamFoam/`)
before doing any of the following:**

- Running a simulation — exact Docker run command and required resource flags (`--shm-size`, `--ulimit`, `--ipc`, `--cpus`, `--memory`)
- Pausing or resuming a run — procedure to signal a clean stop and restore controlDict
- Reconstructing results for ParaView — order of `reconstructParMesh` then `reconstructPar`, and the commands
- Loading results in ParaView — `.foam` marker file, `.pvsm` state file, VTK series file
- Fixing a stale VTK series file after pause/resume — `fix_vtk_series.py` usage
- Deleting processor directories — must use Docker, not `sudo rm`
- Running two cases concurrently — CPU/memory budget guidance
- Making mp4/mpg videos from PNG exports — `ffmpeg` (bundled in the `lbf3` image) collage/concat commands
- Building or rebuilding the Docker image — `CACHE_BUST` pattern
- Post-processing with the synthetic-X-ray script (`results/lateral_xray.py`) or the
  normal-render scripts (`results/lateral_screenshot.py`, `results/top_screenshot.py`) — these run
  in a **different** Docker image (`kitware/paraview:pv-v5.8.0-osmesa-py3`), not `lbf3`; see the
  "Post-processing" section below

Running without the correct Docker flags causes silent failures or poor performance.

Additionally, **always read [vdep_power_sweep.md](vdep_power_sweep.md)** before running, seeding,
or discussing any of `testrun58`–`testrun67` (the active VDEP power-sweep pipeline) — it explains
the seed→fork lineage and which core-count seed each fork must be copied from.

## Repo layout

- `applications/solvers/` — three solvers: `laserbeamFoam`, `compressibleLaserbeamFoam`, `laserMeltFoam`
- `tutorials/laserbeamFoam/vdep/` — VDEP research cases. **Current, in-scope focus is only
  testrun58–67** (see [vdep_power_sweep.md](vdep_power_sweep.md)) — the six bolded forks
  (testrun62–67) are the actual production runs, seeded via testrun58–61.
  testrun30–53 (early prototypes/dead ends) have been removed — see git history for that period
  if needed. testrun54 (a standalone 16→32-core scaling test, not part of the power sweep) and
  testrun55–57 (mesh-resolution tests) still exist on disk but are outside current focus.
- `tutorials/laserbeamFoam/plc/` — PLC reference cases (testrun1–29, 316L steel)
- `tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface/`, `tutorials/multiComponentlaserbeamFoam/` — multi-material cases (`multicomponentLaserbeamFoam` solver)
- `results/lateral_xray.py`, `results/lateral_screenshot.py`, `results/top_screenshot.py`
  — synthetic X-ray and normal-render post-processing scripts for the VDEP power-sweep cases; see
  "Post-processing" below

## Key facts

- Simulations run inside the `lbf3` Docker image (OpenFOAM v2506 + LIGGGHTS + solvers + `ffmpeg`).
  Build with `docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .`
- The correct solver executable is `laserbeamFoam` (not `Flint_multiphaseEulerFoamD`); use `bash -lc` (login shell) inside Docker to source the OpenFOAM profile.
- Processor directories are root-owned — always use Docker (not `sudo rm`) to delete them.
- Always run `reconstructParMesh` before `reconstructPar`.
- After any pause/resume, the VTK series file goes stale. Fix it with: `python3 tutorials/laserbeamFoam/fix_vtk_series.py tutorials/laserbeamFoam/CASE/VTKs`
- Never run a simulation without asking the user first.
- Resuming a paused run needs `log.laserbeamFoam` renamed first (`mv log.laserbeamFoam
  log.laserbeamFoam_stageN`) or `Allrun` silently no-ops — this applies to *any* resume, not just
  multi-stage/core-count-switch runs.

## Post-processing (synthetic X-ray + normal-render views)

- Two **separate** Docker images are involved — don't confuse them:
  - `lbf3` — simulations, `reconstructParMesh`/`reconstructPar`, and `ffmpeg` (video assembly).
  - `kitware/paraview:pv-v5.8.0-osmesa-py3` — headless (OSMesa, no GUI) ParaView/pvpython, used
    only for the `results/*.py` scripts below. Pull with
    `docker pull kitware/paraview:pv-v5.8.0-osmesa-py3` or let the first `docker run` fetch it.
- `lateral_xray.py` — synthetic X-ray attenuation view (lateral/through-thickness),
  with a dotted melt-pool-bottom boundary line derived from a continuous attenuation-ceiling
  threshold (not a boolean presence flag — see the script's own header comment for why that
  matters). `lateral_screenshot.py` — a normal (not ray-traced) ParaView render of the same lateral
  view, colored by lateral position to reveal depth. `top_screenshot.py` — same normal-render
  technique but viewed top-down, colored by height relative to the nominal surface.
  `_render_stacked_video.sh` — drives all three scripts across every reconstructed timestep of a
  given testrun and stacks them into one image per timestep plus an mp4 (see below).
- Invocation (any of the three `results/*.py` scripts takes the same args): `docker run --rm -e
  PYTHONUNBUFFERED=1 -v <repo>:/workspace --entrypoint /opt/paraview/bin/pvpython
  kitware/paraview:pv-v5.8.0-osmesa-py3
  /workspace/results/lateral_xray.py /workspace/<case>.foam <time> <output.png>`.
- Batch-rendering every timestep of a testrun, stacking the three views (top/lateral/X-ray) into one
  image per timestep, and building an mp4: `bash results/_render_stacked_video.sh <testrun>` (bare
  number like `64`, or a full case dir name) -- works for any reconstructed VDEP power-sweep case,
  not just testrun64. Auto-creates a `.foam` marker if the case doesn't have one yet, and errors out
  clearly if the case hasn't been reconstructed (no timestep directories found). Uses the same
  numeric-sorted-timestep + ffmpeg-concat-list pattern as the older single-view
  `results/_render_ls_all_and_video.sh`, since timestep-named files mix scientific/decimal notation
  and don't sort correctly alphabetically.

## Physical parameters (current active cases)

- Material: AlSi10Mg, Tliq=867 K, Tsol=840 K, ρ=2670 kg/m³
- Laser: 1064 nm, 35 µm radius, Drude/Fresnel absorptivity (~11% at normal incidence for liquid Al)
- Process: power varies per fork, 650–900 W effective (see [vdep_power_sweep.md](vdep_power_sweep.md)), 6 m/s scan speed
- Mesh: 40 µm base, AMR at metal–gas interface (3 levels → 5 µm; see e.g.
  `tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al/system/topoSetDict`)
- Decomposition: hierarchical, 32 cores, n=(2 2 8), order yxz
