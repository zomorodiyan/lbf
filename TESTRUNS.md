# Running laserbeamFoam Simulations

All operations run inside the `lbf3` Docker container, which mounts the repo at `/workspace`.

**Case paths:**
- VDEP cases: `tutorials/laserbeamFoam/CASE` (e.g. `testrun30_vdep_1_Al`)
- PLC reference cases: `tutorials/laserbeamFoam/plc/CASE`

Replace `CASE` with the actual case directory name throughout this document.

---

## Prerequisites

### Windows (WSL2)

1. Install **WSL2** if not already present:
   ```powershell
   wsl --install
   ```
2. Install **Docker Desktop** for Windows from [https://www.docker.com/products/docker-desktop/](https://www.docker.com/products/docker-desktop/).
   In Docker Desktop → Settings → Resources → WSL Integration, enable integration for your WSL2 distro.
3. Open a WSL2 terminal and verify: `docker run --rm hello-world`

### Ubuntu / Native Linux

```bash
sudo apt-get update
sudo apt-get install -y docker.io
sudo usermod -aG docker $USER   # log out and back in after this
```

Verify: `docker run --rm hello-world`

---

## Get the repository

Fork the repo on GitHub, then clone your fork:

```bash
git clone https://github.com/YOUR_USERNAME/lbf3.git
cd lbf3
```

---

## Build the Docker image

The solver binary lives inside the image; rebuild after any source changes.
`CACHE_BUST` forces Docker to re-compile while reusing the cached OS and LIGGGHTS layers.

```bash
cd /path/to/lbf3
docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .
```

---

## Run

Most cases use 32 cores:

```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v /path/to/lbf3:/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && bash ./Allrun && echo RUN_COMPLETE"
```

Cases that use 64 cores (e.g. `testrun30_vdep_1_Al`):

```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=64 --memory=76g \
  -v /path/to/lbf3:/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && bash ./Allrun && echo RUN_COMPLETE"
```

Tail the log while it runs:

```bash
tail -f tutorials/laserbeamFoam/CASE/log.laserbeamFoam
```

To run two cases in parallel open two terminals and launch each with a different `CASE`.
Each case uses 32 cores; ensure the host has enough CPUs/RAM for both.

After the run, fix ownership if files are owned by root:

```bash
sudo chown -R "$(id -u):$(id -g)" tutorials/laserbeamFoam/CASE
```

---

## Pause

Signal a clean stop (OpenFOAM finishes the current write before exiting):

```bash
docker exec <container_name> bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && \
   foamDictionary -entry stopAt -set writeNow system/controlDict"
docker wait <container_name>
```

Then restore `endTime` in `system/controlDict`:

```bash
sed -i 's/stopAt.*writeNow;/stopAt          endTime;/' tutorials/laserbeamFoam/CASE/system/controlDict
```

---

## Resume

The Allrun script skips `decomposePar` if `processor0/` already exists, so just rerun the **Run** command above.

---

## Reconstruct for ParaView

**AMR cases: always run `reconstructParMesh` before `reconstructPar`** — it rebuilds the mesh
topology for each timestep. Skipping it causes a FOAM FATAL ERROR on any timestep whose mesh changed.

```bash
docker run --rm -v /path/to/lbf3:/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && \
   reconstructParMesh -noFunctionObjects 2>&1 | tee log.reconstructParMesh && \
   reconstructPar    -noFunctionObjects 2>&1 | tee log.reconstructPar"
```

Both commands skip already-reconstructed timesteps, so re-running after the simulation has advanced further is safe.

---

## Post-processing in ParaView

A ParaView state file is provided at `tutorials/laserbeamFoam/laser.pvsm`.
Open ParaView, go to **File → Load State**, select `laser.pvsm`, and point it at your case directory.
This loads a pre-configured view of the laser, melt pool, and temperature field.

### VTK series files

laserbeamFoam writes a `VTKs/rays_laser0.vtk.series` JSON index alongside each VTK frame.
When a run is **paused and resumed**, OpenFOAM rewrites this file with only the frames from
the new segment — earlier frames disappear from ParaView's timeline.

**Fix:** run the repair script from the repo root:

```bash
# single case
python3 tutorials/laserbeamFoam/fix_vtk_series.py \
  tutorials/laserbeamFoam/CASE/VTKs

# multiple cases at once
python3 tutorials/laserbeamFoam/fix_vtk_series.py \
  tutorials/laserbeamFoam/CASE1/VTKs \
  tutorials/laserbeamFoam/CASE2/VTKs
```

The script scans all `*.vtk` files in the directory, extracts their timestamps from filenames,
and rewrites a complete `.vtk.series` file covering the full simulation history.
You can also ask an AI assistant to do this — point it at the VTKs directory and the script.

---

## Clean up processor directories

`processor*/` directories are large (~25 GB per case). Delete them after a successful reconstruction.
Files are root-owned (written by Docker), so use Docker to remove them:

```bash
docker run --rm -v /path/to/lbf3:/workspace lbf3 bash -c \
  "rm -rf /workspace/tutorials/laserbeamFoam/CASE/processor*/"
```

To clean multiple cases at once:

```bash
docker run --rm -v /path/to/lbf3:/workspace lbf3 bash -c '
for case in testrun27_316L testrun28_316L_PLC; do
  rm -rf /workspace/tutorials/laserbeamFoam/plc/$case/processor*/
  echo "Cleaned $case"
done
'
```

---

## Notes

- Use `cd` inside the container (not the `-w` flag) — the login shell (`bash -lc`) sources the
  OpenFOAM profile which resets `$PWD` to `/root`.
- Processor dirs are root-owned because Docker runs as root; always use Docker (not `sudo rm`) to delete them.
- PLC reference cases live under `tutorials/laserbeamFoam/plc/`; VDEP cases are directly under
  `tutorials/laserbeamFoam/`.
