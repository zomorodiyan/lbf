#!/bin/bash
# Generic version of _render_stacked.sh: instead of 5 hand-picked sample
# timesteps for testrun64, this renders *every* reconstructed timestep of
# a given testrun and assembles the stacked images into an mp4. Works for
# any VDEP power-sweep case (testrun62-67, ...) once it's been reconstructed
# (reconstructParMesh + reconstructPar) and has some timestep directories
# directly under the case folder.
#
# Usage:
#   bash results/_render_stacked_video.sh <testrun>
#     <testrun> is either a bare number (e.g. "64", expanded to
#     testrun64_vdep_3_Al) or a full case directory name
#     (e.g. testrun65_vdep_3_Al).
#
# Output:
#   results/<prefix>_top_screenshot_t<time>.png         (one per timestep)
#   results/<prefix>_lateral_screenshot_t<time>.png     (one per timestep)
#   results/<prefix>_lateral_xray_t<time>.png           (one per timestep)
#   results/<prefix>_transverse_screenshot_t<time>.png  (one per timestep)
#   results/<prefix>_stacked_t<time>.png                (one per timestep)
#   results/<prefix>_stacked_video.mp4
#   where <prefix> is the testrun name, e.g. "testrun65".
set -euo pipefail
cd ~/lbf3

if [ $# -lt 1 ]; then
  echo "Usage: bash results/_render_stacked_video.sh <testrun>  (e.g. 64 or testrun64_vdep_3_Al)"
  exit 1
fi

ARG="$1"
if [[ "$ARG" =~ ^[0-9]+$ ]]; then
  CASE_DIR="testrun${ARG}_vdep_3_Al"
else
  CASE_DIR="$ARG"
fi
PREFIX="${CASE_DIR%%_*}"  # e.g. "testrun65" from "testrun65_vdep_3_Al"

CASE="tutorials/laserbeamFoam/vdep/${CASE_DIR}"
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
echo "Found ${#TIMES[@]} timesteps for $CASE_DIR"

FOAM_FILE=$(find "$CASE" -maxdepth 1 -name '*.foam' | head -1 || true)
if [ -z "$FOAM_FILE" ]; then
  FOAM_FILE="$CASE/case.foam"
  echo "No .foam marker found -- creating one: $FOAM_FILE"
  touch "$FOAM_FILE"
fi
echo "Using .foam marker: $FOAM_FILE"

i=0
for t in "${TIMES[@]}"; do
  i=$((i+1))
  top_png="results/${PREFIX}_top_screenshot_t${t}.png"
  lat_png="results/${PREFIX}_lateral_screenshot_t${t}.png"
  xray_png="results/${PREFIX}_lateral_xray_t${t}.png"
  trans_png="results/${PREFIX}_transverse_screenshot_t${t}.png"
  stacked_png="results/${PREFIX}_stacked_t${t}.png"

  echo "[$i/${#TIMES[@]}] t=$t"

  docker run --rm -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/top_screenshot.py "/workspace/$FOAM_FILE" "$t" "/workspace/$top_png" \
    > /tmp/stackvid_top_${PREFIX}_${t}.log 2>&1 || { echo "  top view FAILED (see /tmp/stackvid_top_${PREFIX}_${t}.log)"; continue; }

  docker run --rm -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/lateral_screenshot.py "/workspace/$FOAM_FILE" "$t" "/workspace/$lat_png" \
    > /tmp/stackvid_lat_${PREFIX}_${t}.log 2>&1 || { echo "  lateral screenshot FAILED (see /tmp/stackvid_lat_${PREFIX}_${t}.log)"; continue; }

  docker run --rm -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/lateral_xray.py "/workspace/$FOAM_FILE" "$t" "/workspace/$xray_png" \
    > /tmp/stackvid_xray_${PREFIX}_${t}.log 2>&1 || { echo "  lateral X-ray FAILED (see /tmp/stackvid_xray_${PREFIX}_${t}.log)"; continue; }

  docker run --rm -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/transverse_screenshot.py "/workspace/$FOAM_FILE" "$t" "/workspace/$trans_png" \
    > /tmp/stackvid_trans_${PREFIX}_${t}.log 2>&1 || { echo "  transverse view FAILED (see /tmp/stackvid_trans_${PREFIX}_${t}.log)"; continue; }

  docker run --rm -v "$(pwd)":/workspace lbf3 bash -lc "
    set -e
    W=\$(for f in /workspace/$top_png /workspace/$lat_png /workspace/$xray_png /workspace/$trans_png; do
      ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 \"\$f\"
    done | sort -n | tail -1)
    ffmpeg -y \
      -i /workspace/$top_png \
      -i /workspace/$lat_png \
      -i /workspace/$xray_png \
      -i /workspace/$trans_png \
      -filter_complex \"[0:v]scale=\${W}:-2[v0];[1:v]scale=\${W}:-2[v1];[2:v]scale=\${W}:-2[v2];[3:v]scale=\${W}:-2[v3];[v0][v1][v2][v3]vstack=inputs=4\" \
      /workspace/$stacked_png -loglevel error
  " > /tmp/stackvid_stack_${PREFIX}_${t}.log 2>&1 || { echo "  stacking FAILED (see /tmp/stackvid_stack_${PREFIX}_${t}.log)"; continue; }
done

echo "ALL FRAMES DONE"

# Build the mp4 in true chronological order via an ffmpeg concat list --
# filenames mix scientific/decimal notation for small early timesteps, so a
# plain alphabetical glob sort would misorder them (same gotcha as
# _render_ls_all_and_video.sh).
CONCAT_LIST="results/_${PREFIX}_stacked_concat.txt"
> "$CONCAT_LIST"
for t in "${TIMES[@]}"; do
  f="results/${PREFIX}_stacked_t${t}.png"
  if [ -f "$f" ]; then
    echo "file '/workspace/$f'" >> "$CONCAT_LIST"
    echo "duration 0.2" >> "$CONCAT_LIST"
  fi
done
# concat demuxer quirk: the last file's duration is ignored, so repeat it once more
LAST_LINE=$(tail -2 "$CONCAT_LIST" | head -1)
echo "$LAST_LINE" >> "$CONCAT_LIST"

VIDEO_OUT="results/${PREFIX}_stacked_video.mp4"
docker run --rm -v "$(pwd)":/workspace lbf3 bash -lc \
  "ffmpeg -y -f concat -safe 0 -i /workspace/$CONCAT_LIST -vsync vfr -pix_fmt yuv420p -c:v libx264 -crf 18 /workspace/$VIDEO_OUT"

echo "VIDEO DONE: $VIDEO_OUT"
