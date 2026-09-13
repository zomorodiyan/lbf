#!/bin/bash
# Driver for measure_melt_vapor_depth.py: resolves a case argument the same
# way as _render_cutaway_updated_video.sh/_render_stacked_video.sh, then runs
# the measurement in the kitware/paraview image (post-processing, not lbf3 --
# see CLAUDE.md's Post-processing section).
#
# Usage:
#   bash results/measure_melt_vapor_depth.sh <case> [extra args passed through
#                                                      to measure_melt_vapor_depth.py]
#     <case> is one of:
#       - a bare number (e.g. "81"), expanded to
#         tutorials/laserbeamFoam/vdep/testrun81_vdep_3_Al
#       - a bare case directory name (e.g. "testrun81_vdep_3_Al")
#       - a path (contains a "/") to any other reconstructed case
#
# Examples:
#   bash results/measure_melt_vapor_depth.sh 81
#   bash results/measure_melt_vapor_depth.sh 84 --steady-frac 0.4 --tmelt solidus
#
# Requires the case to already be reconstructed (reconstructParMesh then
# reconstructPar -- see TESTRUNS.md) before running this.
#
# Output (written to results/, prefix = the case directory's basename):
#   results/<prefix>_melt_vapor_depth.csv
#   results/<prefix>_melt_vapor_depth.png
set -euo pipefail
cd ~/lbf3

if [ $# -lt 1 ]; then
  echo "Usage: bash results/measure_melt_vapor_depth.sh <case> [extra args]"
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
docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace \
  --entrypoint /opt/paraview/bin/pvpython "$IMG" \
  /workspace/results/measure_melt_vapor_depth.py \
  "/workspace/$FOAM_FILE" \
  "/workspace/results/${PREFIX}_melt_vapor_depth.csv" \
  "/workspace/results/${PREFIX}_melt_vapor_depth.png" \
  "$@"

echo "DONE: results/${PREFIX}_melt_vapor_depth.csv, results/${PREFIX}_melt_vapor_depth.png"
