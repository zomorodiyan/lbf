#!/bin/bash
# Generic driver for the 3-panel (top/cutaway/xray) stacked composite --
# sibling to _render_stacked_video.sh's own 4-panel (top/transverse/
# lateral/xray) version, using render_view.py's --view=cutaway (the
# x<-25um half-domain cut, see render_view.py's own header on that view)
# in place of transverse+lateral. Renders every reconstructed timestep of
# a given case and assembles the stacked images into an mp4.
#
# Usage:
#   bash results/_render_cutaway3panel_video.sh <case>
#     <case> is one of:
#       - a bare number (e.g. "69"), expanded to
#         tutorials/laserbeamFoam/vdep/testrun69_vdep_3_Al -- shorthand for the
#         VDEP power-sweep cases only
#       - a bare case directory name (e.g. "testrun69_vdep_3_Al"), also
#         resolved under tutorials/laserbeamFoam/vdep/
#       - a path (contains a "/", relative to repo root or absolute) to any
#         other reconstructed case
#
# Resumable: skips any timestep whose stacked frame already exists, so a
# partial/interrupted run can just be re-invoked. Per-timestep failures
# (either view's render, or the stacking step) are logged and skipped,
# not fatal to the whole batch.
#
# Output (written to results/, prefix = the case directory's basename):
#   results/<prefix>_cutaway3panel_t<time>.png   (one per timestep, stacked)
#   results/<prefix>_cutaway3panel_video.mp4
# Per-view intermediates (top/cutaway/xray PNGs, the cutaway's standalone
# colorbar) are deleted after each timestep's frame is successfully
# stacked -- only the stacked frames and the final video are kept.
set -uo pipefail
cd ~/lbf3

if [ $# -lt 1 ]; then
  echo "Usage: bash results/_render_cutaway3panel_video.sh <case>"
  echo "  <case> is a bare number (69), a bare VDEP case dir name (testrun69_vdep_3_Al),"
  echo "  or a path to any other reconstructed case (tutorials/laserbeamFoam/plc/CASE)."
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
PREFIX="$(basename "$CASE")_cutaway3panel"

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

PARAVIEW_IMG=kitware/paraview:pv-v5.8.0-osmesa-py3

render_view() {
  local view=$1 time=$2 out=$3 log=$4
  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace \
    --entrypoint /opt/paraview/bin/pvpython "$PARAVIEW_IMG" \
    /workspace/results/render_view.py --view="$view" "/workspace/$FOAM_FILE" "$time" "/workspace/$out" \
    > "$log" 2>&1
}

i=0
for t in "${TIMES[@]}"; do
  i=$((i+1))
  stacked_png="results/${PREFIX}_t${t}.png"
  if [ -f "$stacked_png" ]; then
    echo "[$i/${#TIMES[@]}] t=$t already stacked, skipping"
    continue
  fi
  echo "[$i/${#TIMES[@]}] t=$t"

  top_png="results/${PREFIX}_top_t${t}.png"
  cut_png="results/${PREFIX}_cutaway_t${t}.png"
  xray_png="results/${PREFIX}_xray_t${t}.png"
  cut_colorbar="${cut_png%_t*}_colorbar.png"

  render_view top "$t" "$top_png" "/tmp/${PREFIX}_top_${t}.log" \
    || { echo "  top view FAILED (see /tmp/${PREFIX}_top_${t}.log)"; continue; }
  render_view cutaway "$t" "$cut_png" "/tmp/${PREFIX}_cutaway_${t}.log" \
    || { echo "  cutaway view FAILED (see /tmp/${PREFIX}_cutaway_${t}.log)"; continue; }
  render_view xray "$t" "$xray_png" "/tmp/${PREFIX}_xray_${t}.log" \
    || { echo "  xray view FAILED (see /tmp/${PREFIX}_xray_${t}.log)"; continue; }

  docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 \
      bash /workspace/results/_stack_cutaway3panel.sh \
      "$top_png" "$cut_png" "$xray_png" "$stacked_png" "$cut_colorbar" "$t" \
      > "/tmp/${PREFIX}_stack_${t}.log" 2>&1 \
    || { echo "  stacking FAILED (see /tmp/${PREFIX}_stack_${t}.log)"; continue; }

  rm -f "$top_png" "$cut_png" "$xray_png" "$cut_colorbar" \
        "results/${PREFIX}_top_colorbar.png" "results/${PREFIX}_xray_colorbar.png"
  echo "  done: $stacked_png"
done

echo "ALL FRAMES DONE"

CONCAT_LIST="results/_${PREFIX}_concat.txt"
> "$CONCAT_LIST"
LAST_FRAME=""
for t in "${TIMES[@]}"; do
  f="results/${PREFIX}_t${t}.png"
  if [ -f "$f" ]; then
    echo "file '$(basename "$f")'" >> "$CONCAT_LIST"
    echo "duration 0.4" >> "$CONCAT_LIST"
    LAST_FRAME="$f"
  fi
done
if [ -n "$LAST_FRAME" ]; then
  echo "file '$(basename "$LAST_FRAME")'" >> "$CONCAT_LIST"
fi

VIDEO_OUT="results/${PREFIX}_video.mp4"
docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc "
  cd /workspace/results &&
  ffmpeg -y -f concat -safe 0 -i $(basename "$CONCAT_LIST") \
    -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -vsync vfr -pix_fmt yuv420p -c:v libx264 -crf 18 \
    $(basename "$VIDEO_OUT") -loglevel error
"
rm -f "$CONCAT_LIST"

echo "VIDEO DONE: $VIDEO_OUT"
