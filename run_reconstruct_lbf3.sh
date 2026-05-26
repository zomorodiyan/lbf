#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_reconstruct_lbf3.sh
#   ./run_reconstruct_lbf3.sh tutorials/laserbeamFoam/testrun19_316L tutorials/laserbeamFoam/testrun20_316L

CASES=("$@")
if [[ ${#CASES[@]} -eq 0 ]]; then
    CASES=(
        tutorials/laserbeamFoam/testrun19_316L
        tutorials/laserbeamFoam/testrun20_316L
    )
fi

docker run --rm \
  --shm-size=16g \
  --ulimit memlock=-1 \
  --ulimit stack=67108864 \
  --ipc=host \
  --cpus=64 \
  --memory=76g \
  -v "$(pwd):/workspace" \
  -w /workspace \
  lbf3 bash -lc "/workspace/reconstruct_results.sh ${CASES[*]}"
