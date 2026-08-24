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
- Post-processing with `results/render_view.py` (`--view=top|lateral|xray|transverse`) — runs in a
  **different** Docker image (`kitware/paraview:pv-v5.8.0-osmesa-py3`), not `lbf3`; see the
  "Post-processing" section below

Running without the correct Docker flags causes silent failures or poor performance.

Additionally, **always read [vdep_power_sweep.md](vdep_power_sweep.md)** before running, seeding,
or discussing any of `testrun58`–`testrun67` (the active VDEP power-sweep pipeline) — it explains
the seed→fork lineage and which core-count seed each fork must be copied from.

Similarly, **always read [vdep_remedy_sims.md](vdep_remedy_sims.md)** before running, seeding, or
discussing `testrun68`–`testrun70` (the recoil-pressure remedy sims from
[task_Aug4.md](task_Aug4.md)) — it covers the testrun61 seed hand-off reused for these three, the
one-time solver change needed before testrun68, and the Mehrdad/Zixun assignment.

Similarly, **always read [T0_sweep.md](T0_sweep.md)** before running, seeding, or discussing
`testrun73`–`testrun84` (the initial-temperature sweep: four independent testrun60→61→69-style
lineages at T0=100/200/400/500K instead of the usual 300K) — case dictionaries are prepared but
**nothing has been run yet** as of 2026-08-23. Also see
[T0_sweep_assignment.md](T0_sweep_assignment.md) for the Mehrdad/Zixun run split and ordering
(cold pair vs. hot pair) for these cases.

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
- `results/render_view.py` — merged post-processing script (`--view=top|lateral|xray|transverse`)
  for the VDEP power-sweep cases; see "Post-processing" below

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
    only for `results/render_view.py` below. Pull with
    `docker pull kitware/paraview:pv-v5.8.0-osmesa-py3` or let the first `docker run` fetch it.
- `results/render_view.py --view={top,lateral,xray,transverse}` — one merged script, four views
  (used to be 4 separate scripts; merged to kill duplicated code — see the script's own header):
  - `top` — top-down normal render, colored by height relative to the nominal surface.
  - `lateral` — normal render of the lateral (through-thickness) view, colored by lateral position.
  - `xray` — synthetic X-ray attenuation view (same lateral view, different technique: numpy
    ray-tracing + Beer-Lambert, not a ParaView render), with a dotted melt-pool-bottom boundary
    line derived from a continuous attenuation-ceiling threshold (see the script's header for why
    not a boolean presence flag).
  - `transverse` — cross-section perpendicular to the scan direction (looking down the z/scan axis
    at the x-y width-vs-depth plane): three side-by-side panels, each a thin isosurface slab cut at
    a fixed distance (2.0/1.5/1.0mm, left-to-right) behind the laser's current z position.
  - Invocation: `docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace --entrypoint
    /opt/paraview/bin/pvpython kitware/paraview:pv-v5.8.0-osmesa-py3
    /workspace/results/render_view.py --view=<view> /workspace/<case>.foam <time> <output.png>`.
- `results/_render_stacked_video.sh <testrun>` (bare number like `64`, or a full case dir name) —
  batch-renders every reconstructed timestep through all four views, vstacks them into one image
  per timestep, and builds an mp4. Works for any reconstructed VDEP power-sweep case. Auto-creates
  a `.foam` marker if missing; errors out clearly if the case hasn't been reconstructed.

## Physical parameters (current active cases)

- Material: Al 6061, Tliq=925 K, Tsol=855 K, ρ=2700 kg/m³
- Laser: 1064 nm, 35 µm radius, Drude/Fresnel absorptivity (~11% at normal incidence for liquid Al)
- Process: power varies per fork, 650–900 W effective (see [vdep_power_sweep.md](vdep_power_sweep.md)), 6 m/s scan speed
- Mesh: 40 µm base, AMR at metal–gas interface (3 levels → 5 µm; see e.g.
  `tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al/system/topoSetDict`)
- Decomposition: hierarchical, 32 cores, n=(2 2 8), order yxz
