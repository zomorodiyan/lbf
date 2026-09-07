#!/bin/bash
# Generic driver for a single, standalone cutaway-only render across every
# reconstructed timestep of a case, using render_cutaway's current full
# feature set: surface tangential-velocity arrows, cross-section (cut-face)
# velocity arrows, laser ray-tracing overlay, flat-shaded (non-Phong)
# outline tubes, and the bottom-corner-artifact gray fill -- all the
# improvements made this session on top of the plain cutaway view.
# Deliberately kept separate from _render_cutaway3panel_video.sh's own
# cutaway panel for now (user request, 2026-08-27: generate this as its own
# "updated version" first, to evaluate on testrun75 before deciding whether
# to fold these same flags into that 3-panel pipeline's cutaway panel too).
#
# No --highlight: those are hand-tuned (z,y) box coordinates specific to
# testrun69's own melt-pool geometry/timing (see _render_collage_panels.sh)
# and would be meaningless (or wrong) on any other case.
#
# Usage:
#   bash results/_render_cutaway_updated_video.sh <case>
#     <case> is one of:
#       - a bare number (e.g. "75"), expanded to
#         tutorials/laserbeamFoam/vdep/testrun75_vdep_3_Al -- shorthand for
#         the VDEP power-sweep cases only
#       - a bare case directory name (e.g. "testrun75_vdep_3_Al"), also
#         resolved under tutorials/laserbeamFoam/vdep/
#       - a path (contains a "/", relative to repo root or absolute) to any
#         other reconstructed case
#
# Resumable: skips any timestep whose frame PNG already exists.
#
# Output (written to results/, prefix = the case directory's basename):
#   results/<prefix>_cutaway_updated_t<time>.png   (one per timestep)
#   results/<prefix>_cutaway_updated_video.mp4
set -euo pipefail
cd ~/lbf3

if [ $# -lt 1 ]; then
  echo "Usage: bash results/_render_cutaway_updated_video.sh <case>"
  echo "  <case> is a bare number (75), a bare VDEP case dir name (testrun75_vdep_3_Al),"
  echo "  or a path to any other reconstructed case."
  exit 1
fi

ARG="$1"
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

mapfile -t TIMES < <(find "$CASE" -maxdepth 1 -type d -regex ".*/[0-9.e-]+" | xargs -n1 basename | sort -g)
if [ "${#TIMES[@]}" -eq 0 ]; then
  echo "ERROR: no reconstructed timesteps found directly under $CASE"
  echo "  (run reconstructParMesh then reconstructPar first -- see TESTRUNS.md)"
  exit 1
fi
echo "Found ${#TIMES[@]} timesteps for $CASE"

FOAM_FILE=$(find "$CASE" -maxdepth 1 -name '*.foam' | head -1 || true)
if [ -z "$FOAM_FILE" ]; then
  FOAM_FILE="$CASE/case.foam"
  echo "No .foam marker found -- creating one: $FOAM_FILE"
  touch "$FOAM_FILE"
fi
echo "Using .foam marker: $FOAM_FILE"

IMG=kitware/paraview:pv-v5.8.0-osmesa-py3
i=0
for t in "${TIMES[@]}"; do
  i=$((i+1))
  out_png="results/${PREFIX}_cutaway_updated_t${t}.png"
  if [ -f "$out_png" ]; then
    echo "[$i/${#TIMES[@]}] t=$t (already done)"
    continue
  fi
  echo "[$i/${#TIMES[@]}] t=$t"
  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython "$IMG" \
    /workspace/results/render_view.py --view=cutaway \
    --y-min-um=-300 --y-max-um=300 --x-plane-um=-35 --top-crop-frac=0 --supersample=3 --velocity --rays --section-velocity \
    "/workspace/$FOAM_FILE" "$t" "/workspace/$out_png" \
    > /tmp/cutawayupd_${PREFIX}_${t}.log 2>&1 || { echo "  FAILED (see /tmp/cutawayupd_${PREFIX}_${t}.log)"; continue; }
done

echo "ALL FRAMES DONE"

# Build the mp4 in true chronological order via an ffmpeg concat list --
# same scientific/decimal-notation sort gotcha as _render_stacked_video.sh.
CONCAT_LIST="results/_${PREFIX}_cutaway_updated_concat.txt"
> "$CONCAT_LIST"
for t in "${TIMES[@]}"; do
  f="results/${PREFIX}_cutaway_updated_t${t}.png"
  if [ -f "$f" ]; then
    echo "file '/workspace/$f'" >> "$CONCAT_LIST"
    echo "duration 0.8" >> "$CONCAT_LIST"
  fi
done
LAST_LINE=$(tail -2 "$CONCAT_LIST" | head -1)
echo "$LAST_LINE" >> "$CONCAT_LIST"

VIDEO_OUT="results/${PREFIX}_cutaway_updated_video.mp4"
docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc \
  "ffmpeg -y -f concat -safe 0 -i /workspace/$CONCAT_LIST -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -vsync vfr -pix_fmt yuv420p -c:v libx264 -crf 18 /workspace/$VIDEO_OUT"

echo "VIDEO DONE: $VIDEO_OUT"
