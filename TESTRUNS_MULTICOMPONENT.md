# Running multicomponentLaserbeamFoam Simulations

This is the counterpart to [TESTRUNS.md](TESTRUNS.md) for the multi-material solver,
`multicomponentLaserbeamFoam`. For anything not covered here — Docker/Prerequisites setup,
`reconstructParMesh`/`reconstructPar` mechanics, loading results in ParaView, fixing a stale
VTK series file after pause/resume, or making mp4 videos from PNG exports — see TESTRUNS.md;
those procedures are identical for this solver.

All operations run inside the `lbf3` Docker container, which mounts the repo at `/workspace`
(or `/home/mzomoro1/bin/lbf3` for the Cu/LCS cases below — see each case's commands).
All commands below assume you are running from the **repo root** (`cd ~/lbf3` first).

---

## multicomponentLaserbeamFoam — the multi-material solver

`multicomponentLaserbeamFoam` is the solver used for all N-phase dissimilar-material cases
(cases in `tutorials/compressiblelaserbeamFoam/` and `tutorials/multiComponentlaserbeamFoam/`).
It is derived from `compressibleLaserbeamFoam` but treats every gas and vapour phase as
**constant-density** (`rhoConst` EOS, ψ = 0):

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

## Notes

- The active Cu/LCS multi-material case lives under `tutorials/multiComponentlaserbeamFoam/Cu_LCS_interface/` and uses `multicomponentLaserbeamFoam` (3-phase incompressible, no acoustic CFL).
- The older compressible cases under `tutorials/compressiblelaserbeamFoam/` are superseded but kept for reference.
- For general Docker/Prerequisites/ParaView/video-making procedures shared with the single-material solver, see [TESTRUNS.md](TESTRUNS.md).
