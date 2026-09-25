#!/bin/bash
# Driver for measure_melt_length.py: resolves a case argument the same way as
# measure_melt_vapor_depth.sh, then runs the measurement in the kitware/paraview
# image (post-processing, not lbf3 -- see CLAUDE.md's Post-processing section).
#
# Usage:
#   bash results/scripts/measure_melt_length.sh <case> [extra args passed through to
#                                                measure_melt_length.py]
#     <case> is one of:
#       - a bare number (e.g. "81"), expanded to
#         tutorials/laserbeamFoam/vdep/testrun81_vdep_3_Al
#       - a bare case directory name (e.g. "testrun81_vdep_3_Al")
#       - a path (contains a "/") to any other reconstructed case
#
# Examples:
#   bash results/scripts/measure_melt_length.sh 81
#   bash results/scripts/measure_melt_length.sh 84 --steady-start 0.0002 --depth-min-um 10
#   bash results/scripts/measure_melt_length.sh 69 --z-ahead-um 150 --z-behind-um 2000
#
# Requires the case to already be reconstructed (reconstructParMesh then
# reconstructPar -- see TESTRUNS.md) before running this.
#
# Output (written to results/melt_pool_measurements/, prefix = the case directory's basename):
#   results/melt_pool_measurements/<prefix>_melt_length.csv
#   results/melt_pool_measurements/<prefix>_melt_length.png
set -euo pipefail
# repo root = parent of this script's dir, so it works on any checkout
# (was hardcoded to one machine's path)
cd "$(dirname "$(readlink -f "$0")")/../.."

if [ $# -lt 1 ]; then
  echo "Usage: bash results/scripts/measure_melt_length.sh <case> [extra args]"
  echo "  <case> is a bare number (81), a bare VDEP case dir name (testrun81_vdep_3_Al),"
  echo "  or a path to any other reconstructed case."
  exit 1
fi

ARG="$1"
shift
if [[ "$ARG" =~ ^[0-9]+$ ]]; then
  CASE="tutorials/laserbeamFoam/vdep/testrun${ARG}_vdep_3_Al"
elif [[ "$ARG" == */* ]]; then
  CASE="${ARG%/}"
else
  CASE="tutorials/laserbeamFoam/vdep/${ARG}"
fi
PREFIX="$(basename "$CASE")"

if [ ! -d "$CASE" ]; then
  echo "ERROR: case directory not found: $CASE"
  exit 1
fi

FOAM_FILE=$(find "$CASE" -maxdepth 1 -name '*.foam' | head -1 || true)
if [ -z "$FOAM_FILE" ]; then
  FOAM_FILE="$CASE/case.foam"
  echo "No .foam marker found -- creating one: $FOAM_FILE"
  touch "$FOAM_FILE"
fi
echo "Using .foam marker: $FOAM_FILE"

IMG=kitware/paraview:pv-v5.8.0-osmesa-py3
mkdir -p results/melt_pool_measurements
docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace \
  --entrypoint /opt/paraview/bin/pvpython "$IMG" \
  /workspace/results/scripts/measure_melt_length.py \
  "/workspace/$FOAM_FILE" \
  "/workspace/results/melt_pool_measurements/${PREFIX}_melt_length.csv" \
  "/workspace/results/melt_pool_measurements/${PREFIX}_melt_length.png" \
  "$@"

echo "DONE: results/melt_pool_measurements/${PREFIX}_melt_length.csv, results/melt_pool_measurements/${PREFIX}_melt_length.png"
