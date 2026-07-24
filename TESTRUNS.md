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

## Multi-material solver (multicomponentLaserbeamFoam)

Run/pause/reconstruct procedures for `multicomponentLaserbeamFoam` (the N-phase dissimilar-
material solver used by `SS316L_Ti64_interface` and the `Cu_LCS_*` cases) have moved to
[TESTRUNS_MULTICOMPONENT.md](TESTRUNS_MULTICOMPONENT.md) — read that file instead of this one
before running, pausing, or reconstructing any case that uses that solver.

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
- Multi-material cases (`tutorials/multiComponentlaserbeamFoam/`, `tutorials/compressiblelaserbeamFoam/`) use the `multicomponentLaserbeamFoam` solver — see [TESTRUNS_MULTICOMPONENT.md](TESTRUNS_MULTICOMPONENT.md) for their run procedures.
