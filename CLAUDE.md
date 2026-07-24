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
- Making mp4 videos from PNG exports — `ffmpeg` collage commands
- Building or rebuilding the Docker image — `CACHE_BUST` pattern

Running without the correct Docker flags causes silent failures or poor performance.

Additionally, **always read [vdep_power_sweep.md](vdep_power_sweep.md)** before running, seeding,
or discussing any of `testrun58`–`testrun67` (the active VDEP power-sweep pipeline) — it explains
the seed→fork lineage and which core-count seed each fork must be copied from.

## Repo layout

- `applications/solvers/` — three solvers: `laserbeamFoam`, `compressibleLaserbeamFoam`, `laserMeltFoam`
- `tutorials/laserbeamFoam/vdep/` — active VDEP research cases (testrun30+). Current production
  pipeline is the **testrun58–67 power sweep** (see [vdep_power_sweep.md](vdep_power_sweep.md));
  testrun35/36 (AlSi10Mg, 700 W, 6 m/s) is the earlier single-case baseline, still referenced for
  physical parameters below.
- `tutorials/laserbeamFoam/plc/` — PLC reference cases (testrun1–29, 316L steel)
- `tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface/`, `tutorials/multiComponentlaserbeamFoam/` — multi-material cases (`multicomponentLaserbeamFoam` solver)

## Key facts

- Everything runs inside the `lbf3` Docker image (OpenFOAM v2506 + LIGGGHTS + solvers). Build with `docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .`
- The correct solver executable is `laserbeamFoam` (not `Flint_multiphaseEulerFoamD`); use `bash -lc` (login shell) inside Docker to source the OpenFOAM profile.
- Processor directories are root-owned — always use Docker (not `sudo rm`) to delete them.
- Always run `reconstructParMesh` before `reconstructPar`.
- After any pause/resume, the VTK series file goes stale. Fix it with: `python3 tutorials/laserbeamFoam/fix_vtk_series.py tutorials/laserbeamFoam/CASE/VTKs`
- Never run a simulation without asking the user first.

## Physical parameters (current active cases)

- Material: AlSi10Mg, Tliq=867 K, Tsol=840 K, ρ=2670 kg/m³
- Laser: 1064 nm, 35 µm radius, Drude/Fresnel absorptivity (~11% at normal incidence for liquid Al)
- Process: 700 W effective (1000 W nominal × 70%), 6 m/s scan speed
- Mesh: 40 µm base, AMR at metal–gas interface (testrun35: 2 levels → 10 µm; testrun36: 3 levels → 5 µm)
- Decomposition: hierarchical, 32 cores, n=(2 2 8), order yxz
