#!/bin/bash
set -e
cd /workspace/tutorials/laserbeamFoam/vdep/testrun55_vdep_mesh_test
rm -rf constant/polyMesh log.* 0 constant/cellSets constant/refineHexMesh* 2>/dev/null || true

echo "=== blockMesh ==="
time blockMesh > log.blockMesh 2>&1
tail -5 log.blockMesh

echo "=== checkMesh (base) ==="
checkMesh > log.checkMesh.base 2>&1
grep -E "cells:|Mesh OK|Failed" log.checkMesh.base

for i in 1 2 3; do
  echo "=== refine pass $i ==="
  time topoSet -dict system/topoSetDict > log.topoSet.$i 2>&1
  grep -E "cells|Set" log.topoSet.$i | tail -3
  time refineHexMesh boxToRefine -overwrite > log.refineHexMesh.$i 2>&1
  tail -5 log.refineHexMesh.$i
done

echo "=== checkMesh (final) ==="
checkMesh > log.checkMesh.final 2>&1
grep -E "cells:|Mesh OK|Failed" log.checkMesh.final

echo "=== cgroup peak memory ==="
cat /sys/fs/cgroup/memory.peak
