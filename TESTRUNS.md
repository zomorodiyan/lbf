# Running and Reconstructing Testruns

All operations run inside the `lbf3` Docker container, which mounts the repo at `/workspace`.
Replace `CASE` with the actual case directory name (e.g. `testrun27_316L`).

---

## Rebuild Docker image (required after source changes)

The solver binary lives inside the Docker image; source changes require a full image rebuild.
The `CACHE_BUST` argument forces Docker to re-run the `COPY` + compile steps while reusing
the cached OS and LIGGGHTS layers (much faster than `--no-cache`).

```bash
cd /home/mzomoro1/bin/lbf3
docker build --build-arg CACHE_BUST=$(date +%s) -t lbf3 .
```

---

## Run

```bash
docker run --rm --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=76g \
  -v /home/mzomoro1/bin/lbf3:/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && bash ./Allrun && echo RUN_COMPLETE"
```

Tail the log while it runs:

```bash
tail -f tutorials/laserbeamFoam/CASE/log.laserbeamFoam
```

To run two cases in parallel, open two terminals and launch each with different `CASE` names.
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
docker run --rm -v /home/mzomoro1/bin/lbf3:/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/CASE && \
   reconstructParMesh -noFunctionObjects 2>&1 | tee log.reconstructParMesh && \
   reconstructPar    -noFunctionObjects 2>&1 | tee log.reconstructPar"
```

Both commands skip already-reconstructed timesteps, so it is safe to re-run after the simulation has advanced further.

---

## Clean up processor directories

`processor*/` directories are large (~25G per case). Delete them after a successful reconstruction.
Files are owned by root (written by Docker), so use Docker to remove them:

```bash
docker run --rm -v /home/mzomoro1/bin/lbf3:/workspace lbf3 bash -c \
  "rm -rf /workspace/tutorials/laserbeamFoam/CASE/processor*/"
```

To clean multiple cases at once:

```bash
docker run --rm -v /home/mzomoro1/bin/lbf3:/workspace lbf3 bash -c '
for case in testrun27_316L testrun28_316L_PLC; do
  rm -rf /workspace/tutorials/laserbeamFoam/$case/processor*/
  echo "Cleaned $case"
done
'
```

---

## Notes

- Use `cd` inside the container (not the `-w` flag) — the login shell (`bash -lc`) sources the OpenFOAM profile which resets `$PWD` to `/root`.
- Processor dirs are root-owned because Docker runs as root; always use Docker (not `sudo rm`) to delete them.
