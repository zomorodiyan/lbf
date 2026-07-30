#!/bin/bash
set -euo pipefail
cd /home/mzomoro1/lbf3

docker run -d --name testrun66_run --shm-size=32g --ulimit memlock=-1 --ulimit stack=67108864 \
  --ipc=host --cpus=32 --memory=56g \
  -v "$(pwd)":/workspace lbf3 bash -lc \
  "cd /workspace/tutorials/laserbeamFoam/vdep/testrun66_vdep_3_Al && bash ./Allrun && echo RUN_COMPLETE"

echo "container started:"
docker ps -a --filter name=testrun66_run
