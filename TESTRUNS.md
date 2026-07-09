# Running laserbeamFoam Simulations

All operations run inside the `lbf3` Docker container, which mounts the repo at `/workspace`.
All commands below assume you are running from the **repo root** (`cd ~/lbf3` first).

**Case paths:**
- VDEP cases: `tutorials/laserbeamFoam/vdep/CASE` — e.g. `testrun32_vdep_3_Al`
- PLC reference cases: `tutorials/laserbeamFoam/plc/CASE`

Replace `CASE` with the actual case directory name throughout this document.

---

## Prerequisites

### System requirements

- Ubuntu 20.04 or 22.04
- At least **32 CPU cores** and **64 GB RAM** (cases are configured for 32 cores and ~76 GB)
- ~50 GB free disk space per case (processor directories are large)

### Docker

```bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo usermod -aG docker $USER
newgrp docker          # apply group change without logging out
docker run --rm hello-world   # verify
```

### ParaView (for post-processing)

```bash
sudo apt-get install -y paraview
```

---

## Get the repository

Fork [Mehrdad's repo](https://github.com/zomorodiyan/lbf) on GitHub, then clone your fork:

```bash
git clone https://github.com/YOUR_USERNAME/lbf.git ~/lbf3
cd ~/lbf3
```

---

## How the Docker setup works

Rather than installing OpenFOAM and the custom solver directly on your machine, everything runs
inside a **Docker container** built from the `Dockerfile` at the repo root.

The image contains:
- OpenFOAM v2506 (the CFD framework)
- LIGGGHTS (discrete element method solver, used by powder-bed cases)
- The compiled `laserbeamFoam` solver binary (built from `applications/` in this repo)

When you run a simulation, Docker mounts the repo directory into the container at `/workspace`
so the case files on your machine are directly accessible — results are written back to your
disk and survive after the container exits.

## Build the Docker image

Build once after cloning (and again after any source code change in `applications/`).
`CACHE_BUST` forces Docker to re-run the compile step while reusing the cached OS and LIGGGHTS
layers, making rebuilds much faster than a full `--no-cache` build.

```bash
cd ~/lbf3
docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .
```

This takes ~10–20 minutes the first time, ~2–5 minutes for subsequent rebuilds.

---

## Recommended starting case

**`vdep/testrun35_vdep_3_Al`** — AlSi10Mg, 700 W effective, 6 m/s, 35 µm beam radius, 32 cores, 2-level AMR (→10 µm).
**`vdep/testrun36_vdep_3_Al`** — same parameters, 3-level AMR (→5 µm), 25% longer track. Current active case.

---

## Run

```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && bash ./Allrun && echo RUN_COMPLETE"
```

Adjust `--memory=76g` down if your workstation has less RAM (minimum ~48g for 32-core cases).

Tail the log while it runs (in a second terminal):

```bash
tail -f tutorials/laserbeamFoam/CASE/log.laserbeamFoam
```

After the run, fix ownership if files are owned by root:

```bash
sudo chown -R "$(id -u):$(id -g)" tutorials/laserbeamFoam/CASE
```

---

## Pause

Signal a clean stop (OpenFOAM finishes the current timestep before exiting):

```bash
docker exec <container_name> bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && \
   foamDictionary -entry stopAt -set writeNow system/controlDict"
docker wait <container_name>
```

Find `<container_name>` with `docker ps`.

Then restore `endTime` so the case can be resumed:

```bash
sed -i 's/stopAt.*writeNow;/stopAt          endTime;/' \
  tutorials/laserbeamFoam/CASE/system/controlDict
```

---

## Resume

The Allrun script skips `decomposePar` if `processor0/` already exists, so just rerun the **Run** command above.

**Multi-stage runs** (e.g. testrun36, testrun43 — different write intervals or power profiles per stage):

1. Set `endTime` and `writeInterval` in `system/controlDict` for stage 1 and run normally.
2. When stage 1 finishes, **rename the log file** to preserve it:
   ```bash
   mv tutorials/laserbeamFoam/CASE/log.laserbeamFoam \
      tutorials/laserbeamFoam/CASE/log.laserbeamFoam_stage1
   ```
3. Update `endTime` and `writeInterval` (and any other parameters, e.g. `timeVsLaserPower`) for stage 2.
4. Rerun the same **Run** command — OpenFOAM restarts from `latestTime` automatically, `decomposePar` is skipped because `processor0/` already exists.
5. Repeat for further stages, incrementing the suffix (`_stage2`, `_stage3`, …).

**Example — testrun43** (dense snapshots during violent startup, coarser after):
- Stage 1: `endTime 1e-5`, `writeInterval 1e-6` (0–10 µs, every 1 µs)
- Stage 2: `endTime 1e-4`, `writeInterval 4e-6` (10–100 µs, every 4 µs)

---

## Reconstruct for ParaView

**Always run `reconstructParMesh` before `reconstructPar`** — it rebuilds the mesh topology
for each timestep. Skipping it causes a FOAM FATAL ERROR on any timestep whose mesh changed.

```bash
docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && \
   reconstructParMesh -noFunctionObjects 2>&1 | tee log.reconstructParMesh && \
   reconstructPar    -noFunctionObjects 2>&1 | tee log.reconstructPar"
```

Both commands skip already-reconstructed timesteps, so re-running after the simulation has
advanced further is safe.

---

## Post-processing in ParaView

A pre-configured ParaView state is at `tutorials/laserbeamFoam/laser.pvsm`.
It stores the camera, filters, and colour maps — you just redirect it to your case files.

### Step-by-step: loading the state for a new case

1. Make sure the case has been reconstructed (see above) and a `.foam` marker file exists.
   If it doesn't, create one:
   ```bash
   touch tutorials/laserbeamFoam/CASE/CASE.foam
   ```

2. Open ParaView → **File → Load State** → navigate to `tutorials/laserbeamFoam/laser.pvsm` → **OK**.

3. In the dialog that appears, select **"Choose File Names"** (not "Use File Names From State" —
   that points to the original machine's paths and will fail).

4. ParaView lists each data source that needs a file. Point them to your case:
   - **OpenFOAM reader** → browse to `tutorials/laserbeamFoam/CASE/CASE.foam`
   - **VTK series** → browse to `tutorials/laserbeamFoam/CASE/VTKs/rays_laser0.vtk.series`

5. Click **OK**. The melt pool, temperature field, and laser rays load automatically.

### VTK series files — fixing after pause/resume

laserbeamFoam writes a `VTKs/rays_laser0.vtk.series` JSON index alongside each VTK frame.
When a run is **killed or paused and resumed**, OpenFOAM rewrites this file with only the frames
from the new segment — earlier frames disappear from ParaView's timeline.

**Fix:** run the repair script from the repo root:

```bash
python3 tutorials/laserbeamFoam/fix_vtk_series.py tutorials/laserbeamFoam/CASE/VTKs
```

Alternatively, ask an AI assistant — point it at `fix_vtk_series.py` and the `VTKs/` directory
and it will repair the file automatically.

---

## Clean up processor directories

`processor*/` directories are large (~25 GB per case). Delete them after a successful reconstruction.
Files are root-owned (written by Docker), so always use Docker to remove them:

```bash
docker run --rm -v $(pwd):/workspace lbf3 bash -c \
  "rm -rf /workspace/tutorials/laserbeamFoam/CASE/processor*/"
```

---

## multicomponentLaserbeamFoam — the multi-material solver

`multicomponentLaserbeamFoam` is the solver used for all N-phase dissimilar-material cases
(cases in `tutorials/compressiblelaserbeamFoam/`). It is derived from `compressibleLaserbeamFoam`
but treats every gas and vapour phase as **constant-density** (`rhoConst` EOS, ψ = 0):

| Feature | `laserbeamFoam` | `multicomponentLaserbeamFoam` |
|---|---|---|
| Number of metal phases | 1 | N (SS316L + Ti64, Cu + steel, …) |
| Gas EOS | constant ρ | constant ρ (rhoConst) |
| Pressure equation | incompressible Poisson | incompressible Poisson |
| Acoustic CFL constraint | none | none |
| Compressible work in TEqn | no | no |

### Build (required after any source change)

The Dockerfile builds all solvers automatically via `./Allwmake`. After modifying source code,
rebuild the image:

```bash
cd ~/lbf3
docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .
```

To build only this solver inside a running container (faster for iteration):

```bash
docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/applications/solvers/multicomponentLaserbeamFoam && ./Allwmake"
```

---

## Multi-material case: SS316L / Ti-6Al-4V interface

Case path: `tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface`

Uses `multicomponentLaserbeamFoam` (5-phase: SS316L + vapour, Ti64 + vapour, air).
Domain: 400 × 160 × 200 µm. Base mesh: 40 µm. 3 AMR levels → 5 µm finest.
Laser: 200 W, 1000 mm/s, 35 µm radius, scans +x across the 316L/Ti64 interface.
Physical time: 300 µs. Decomposed across 32 cores (hierarchical, n=2×2×8).

### Run

```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface && bash ./Allrun && echo RUN_COMPLETE"
```

Tail the log:

```bash
tail -f tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface/log.multicomponentLaserbeamFoam
```

### Pause

```bash
docker exec <container_name> bash -lc \
  "cd /workspace/tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface && \
   foamDictionary -entry stopAt -set writeNow system/controlDict"
docker wait <container_name>
```

Restore after pause:

```bash
sed -i 's/stopAt.*writeNow;/stopAt          endTime;/' \
  tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface/system/controlDict
```

### Reconstruct

```bash
docker run --rm -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface && \
   reconstructParMesh -noFunctionObjects 2>&1 | tee log.reconstructParMesh && \
   reconstructPar    -noFunctionObjects 2>&1 | tee log.reconstructPar"
```

### Running concurrently with another case

The two containers are independent processes — Docker imposes the `--cpus` and `--memory` limits
per container, so as long as the combined totals fit your hardware, both run simultaneously.
Example: two 32-core cases on a 64-core machine:

```bash
# Terminal 1 — existing vdep run (already running, or launch it):
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/vdep/CASE && bash ./Allrun && echo RUN_COMPLETE"

# Terminal 2 — new multi-material case:
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v $(pwd):/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface && bash ./Allrun && echo RUN_COMPLETE"
```

Both write to separate directories under `/workspace` so there is no file conflict.

---

## Multi-material case: Cu / Low-Carbon Steel (Soysal & Kou replication)

Case path: `tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface`

Uses `multicomponentLaserbeamFoam` (3-phase: Cu + steel + air).
Replicates Soysal & Kou 2016 (Acta Materialia) — macrosegregation in Cu/steel GTAW butt welding.
See [tasks_multi.md](tasks_multi.md) for full physics context and progress tracking.

**Domain:** x=5 mm (Cu: 0–4 mm, steel: 4–5 mm), y=2 mm (gas/metal/gas), z=6 mm (scan)
**Mesh:** 200 µm base, no AMR. 25×10×30 = 7500 base cells.
**Laser:** 200 W nominal (test) / 4600 W (production), 500 µm radius, conduction mode.
**Beam center:** x=2.5 mm (1.5 mm offset into Cu from interface at x=4 mm), y=1.7 mm, scans +z at 1 mm/s.
**Physical time:** 3 s (z=1.5 mm to z=4.5 mm). maxDeltaT=1e-4 s.
**Decomposition:** 8 cores, hierarchical n=(4 1 2).
**Timestep:** ~1e-4 s at steady state — 180× faster than compressibleLaserbeamFoam (no acoustic CFL).

### Run (parallel, 8 cores)

```bash
docker run -d --name cu_lcs \
  --shm-size=8g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=8 --memory=20g \
  -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  -v $(pwd):/home/mzomoro1/bin/lbf3 lbf3 bash -lc \
  "cd /home/mzomoro1/bin/lbf3/tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface && \
   rm -rf processor* log.* 0 && cp -r initial 0 && \
   blockMesh > log.blockMesh 2>&1 && \
   setFields > log.setFields 2>&1 && \
   decomposePar > log.decomposePar 2>&1 && \
   mpirun --allow-run-as-root -np 8 multicomponentLaserbeamFoam -parallel \
   > log.multicomponentLaserbeamFoam 2>&1 && echo RUN_COMPLETE"
```

Tail the log:

```bash
tail -f tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/log.multicomponentLaserbeamFoam
```

### Pause

```bash
docker exec cu_lcs bash -lc \
  "foamDictionary -entry stopAt -set writeNow \
   /home/mzomoro1/bin/lbf3/tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/system/controlDict"
docker wait cu_lcs
```

Restore after pause:

```bash
sed -i 's/stopAt.*writeNow;/stopAt          endTime;/' \
  tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/system/controlDict
```

### Reconstruct

```bash
docker run --rm \
  -e OMPI_ALLOW_RUN_AS_ROOT=1 -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  -v $(pwd):/home/mzomoro1/bin/lbf3 lbf3 bash -lc \
  "cd /home/mzomoro1/bin/lbf3/tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface && \
   reconstructParMesh -noFunctionObjects 2>&1 | tee log.reconstructParMesh && \
   reconstructPar    -noFunctionObjects 2>&1 | tee log.reconstructPar"
```

**Notes on reconstruction output:**
- `reconstructParMesh` prints `No meshes` — this is normal for a static (non-AMR) mesh; the base mesh in `constant/polyMesh` is already reconstructed and no per-timestep mesh is written.
- `reconstructPar` prints `No FA fields` — normal; no face-agglomeration fields in this case.
- `DilationError` is a **solver-written field** (volumetric dilation from phase change), not an error message. It will appear in every time directory and in `log.reconstructPar` output — ignore it.
- After reconstruction, all time directories (0.05, 0.1, ...) appear directly under the case root alongside `processor*/`.
- Open in ParaView using the `.foam` marker file: `Cu_LCS_interface/Cu_LCS.foam`
- Key fields to visualise: `alpha.Cu`, `alpha.steel`, `T`, `U`, `LatentHeat` (non-zero where metal is molten)

### Clean processor dirs

```bash
docker run --rm -v $(pwd):/home/mzomoro1/bin/lbf3 lbf3 bash -c \
  "rm -rf /home/mzomoro1/bin/lbf3/tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/processor*/"
```

---

## Making videos from ParaView PNG exports

ParaView exports image sequences as `name.0000.png`, `name.0001.png`, etc. Use `ffmpeg` to combine them into an mp4.

**Single view (e.g. top view only):**
```bash
ffmpeg -y -framerate 5 -i tutorials/laserbeamFoam/vdep/CASE/top.%04d.png \
  -c:v libx264 -pix_fmt yuv420p -crf 18 \
  tutorials/laserbeamFoam/vdep/CASE/out.mp4
```

**Two views stacked vertically (top above bottom):**
```bash
ffmpeg -y -framerate 5 -i tutorials/laserbeamFoam/vdep/CASE/top.%04d.png \
  -framerate 5 -i tutorials/laserbeamFoam/vdep/CASE/bot.%04d.png \
  -filter_complex "[0:v][1:v]vstack=inputs=2" \
  -c:v libx264 -pix_fmt yuv420p -crf 18 -r 5 \
  tutorials/laserbeamFoam/vdep/CASE/out.mp4
```

**Two views side by side:**
```bash
ffmpeg -y -framerate 5 -i tutorials/laserbeamFoam/vdep/CASE/left.%04d.png \
  -framerate 5 -i tutorials/laserbeamFoam/vdep/CASE/right.%04d.png \
  -filter_complex "[0:v][1:v]hstack=inputs=2" \
  -c:v libx264 -pix_fmt yuv420p -crf 18 -r 5 \
  tutorials/laserbeamFoam/vdep/CASE/out.mp4
```

Notes:
- Adjust `--framerate` to control playback speed (5 fps is good for melt pool animations).
- `-crf 18` gives high quality; increase to 23–28 for smaller files.
- If frames have different resolutions, add `-vf scale=WIDTH:HEIGHT` before the output to normalise.
- Output is at `\\wsl$\Ubuntu\home\mzomoro1\bin\lbf3\...` when browsing from Windows.

---

## Notes

- Run all `docker run` commands from the repo root (`~/lbf3`) so `$(pwd)` resolves correctly.
- Use `cd` inside the container (not the `-w` flag) — the login shell (`bash -lc`) sources the
  OpenFOAM profile which resets `$PWD` to `/root`.
- Processor dirs are root-owned because Docker runs as root; always use Docker (not `sudo rm`) to delete them.
- PLC reference cases (testruns 1–29) live under `tutorials/laserbeamFoam/plc/`.
- VDEP cases (testruns 30+) live under `tutorials/laserbeamFoam/vdep/`.
- The active Cu/LCS multi-material case lives under `tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/` and uses `multicomponentLaserbeamFoam` (3-phase incompressible, no acoustic CFL).
- The older compressible cases under `tutorials/compressiblelaserbeamFoam/` are superseded but kept for reference.
