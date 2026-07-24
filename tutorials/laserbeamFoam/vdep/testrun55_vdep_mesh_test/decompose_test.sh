#!/bin/bash
set -e
cd /workspace/tutorials/laserbeamFoam/vdep/testrun55_vdep_mesh_test
rm -rf processor* 2>/dev/null || true
echo "=== decomposePar ==="
time decomposePar > log.decomposePar 2>&1
tail -8 log.decomposePar
echo "=== cgroup peak memory ==="
cat /sys/fs/cgroup/memory.peak
