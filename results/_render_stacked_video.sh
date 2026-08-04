#!/bin/bash
# Generic version of _render_stacked.sh: instead of 5 hand-picked sample
# timesteps for testrun64, this renders *every* reconstructed timestep of
# a given case and assembles the stacked images into an mp4. Works for any
# reconstructed OpenFOAM case (not just the VDEP power-sweep cases) once it's
# been reconstructed (reconstructParMesh + reconstructPar) and has some
# timestep directories directly under the case folder.
#
# Usage:
#   bash results/_render_stacked_video.sh <case>
#     <case> is one of:
#       - a bare number (e.g. "64"), expanded to
#         tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al -- shorthand for the
#         VDEP power-sweep cases only
#       - a bare case directory name (e.g. "testrun65_vdep_3_Al"), also
#         resolved under tutorials/laserbeamFoam/vdep/
#       - a path (contains a "/", relative to repo root or absolute) to any
#         other reconstructed case, e.g.
#         tutorials/laserbeamFoam/plc/testrun12_1_SS316L or
#         /workspace/tutorials/compressiblelaserbeamFoam/SS316L_Ti64_interface
#
# Drives results/render_view.py (--view=top/lateral/xray/transverse) --
# formerly 4 separate scripts, merged into one (see that file's own header).
# Three passes: render every view for every timestep, then normalize each
# view's own time series to one shared size (render_view.py's
# --trim-batch, see its own docstring) since each render_* function's
# per-frame content trim otherwise leaves that view's dimensions drifting
# timestep to timestep (most visibly render_lateral, whose content height
# tracks the melt-pool/spatter depth -- user-reported, 2026-08-03), then
# stack.
#
# Output (written to results/, prefix = the case directory's basename, e.g.
# "testrun65_vdep_3_Al" or "testrun12_1_SS316L" -- not shortened, so cases
# from different directories never collide). Output filenames keep their
# original per-view naming (top_screenshot/lateral_screenshot/lateral_xray/
# transverse_screenshot) even though render_view.py's own --view flags are
# shorter (top/lateral/xray/transverse) -- unrelated concerns, this script
# controls the output filenames itself, independent of how the render
# script names its own view modes:
#   results/<prefix>_top_screenshot_t<time>.png         (one per timestep --
#                                                          colorbar embedded
#                                                          inside, see
#                                                          render_view.py's
#                                                          _save_colorbar())
#   results/<prefix>_lateral_screenshot_t<time>.png     (one per timestep --
#                                                          colorbar embedded
#                                                          inside, ditto)
#   results/<prefix>_lateral_xray_t<time>.png           (one per timestep)
#   results/<prefix>_transverse_screenshot_t<time>.png  (one per timestep)
#   results/<prefix>_stacked_t<time>.png                (one per timestep --
#                                                          top/transverse/
#                                                          lateral/xray
#                                                          vstacked, in that
#                                                          order)
#   results/<prefix>_stacked_video.mp4
set -euo pipefail
cd ~/lbf3

if [ $# -lt 1 ]; then
  echo "Usage: bash results/_render_stacked_video.sh <case>"
  echo "  <case> is a bare number (64), a bare VDEP case dir name (testrun64_vdep_3_Al),"
  echo "  or a path to any other reconstructed case (tutorials/laserbeamFoam/plc/CASE)."
  exit 1
fi

ARG="$1"
if [[ "$ARG" =~ ^[0-9]+$ ]]; then
  CASE="tutorials/laserbeamFoam/vdep/testrun${ARG}_vdep_3_Al"
elif [[ "$ARG" == */* ]]; then
  CASE="${ARG%/}"  # strip a trailing slash, if any
else
  CASE="tutorials/laserbeamFoam/vdep/${ARG}"
fi
PREFIX="$(basename "$CASE")"  # e.g. "testrun65_vdep_3_Al" or "testrun12_1_SS316L"

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

# Pass 1: render every view for every timestep first (no stacking yet).
# Each render_* function trims its own PNG down to *that frame's own*
# content extent (_trim_side_whitespace/_trim_vertical_whitespace) -- fine
# within a single frame, but the content extent (especially render_lateral's
# melt-pool/spatter depth) genuinely differs timestep to timestep, so each
# view's per-frame dimensions drifted across a batch and so did the final
# stacked frame size (user-reported, 2026-08-03). Pass 2 below normalizes
# each view's own time series to one shared size before anything gets
# stacked, so collect this pass's per-view file lists as we go.
top_pngs=(); lat_pngs=(); xray_pngs=(); trans_pngs=()
i=0
for t in "${TIMES[@]}"; do
  i=$((i+1))
  top_png="results/${PREFIX}_top_screenshot_t${t}.png"
  lat_png="results/${PREFIX}_lateral_screenshot_t${t}.png"
  xray_png="results/${PREFIX}_lateral_xray_t${t}.png"
  trans_png="results/${PREFIX}_transverse_screenshot_t${t}.png"

  echo "[$i/${#TIMES[@]}] t=$t"

  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/render_view.py --view=top "/workspace/$FOAM_FILE" "$t" "/workspace/$top_png" \
    > /tmp/stackvid_top_${PREFIX}_${t}.log 2>&1 || { echo "  top view FAILED (see /tmp/stackvid_top_${PREFIX}_${t}.log)"; continue; }

  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/render_view.py --view=lateral "/workspace/$FOAM_FILE" "$t" "/workspace/$lat_png" \
    > /tmp/stackvid_lat_${PREFIX}_${t}.log 2>&1 || { echo "  lateral screenshot FAILED (see /tmp/stackvid_lat_${PREFIX}_${t}.log)"; continue; }

  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/render_view.py --view=xray "/workspace/$FOAM_FILE" "$t" "/workspace/$xray_png" \
    > /tmp/stackvid_xray_${PREFIX}_${t}.log 2>&1 || { echo "  lateral X-ray FAILED (see /tmp/stackvid_xray_${PREFIX}_${t}.log)"; continue; }

  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/render_view.py --view=transverse "/workspace/$FOAM_FILE" "$t" "/workspace/$trans_png" \
    > /tmp/stackvid_trans_${PREFIX}_${t}.log 2>&1 || { echo "  transverse view FAILED (see /tmp/stackvid_trans_${PREFIX}_${t}.log)"; continue; }

  top_pngs+=("$top_png"); lat_pngs+=("$lat_png"); xray_pngs+=("$xray_png"); trans_pngs+=("$trans_png")
done

echo "ALL FRAMES RENDERED"

# Pass 2: normalize each view's own time series to one shared size by
# padding every frame up to that view's own max width/height across the
# whole batch (see render_view.py's _pad_to_uniform_size docstring for why
# padding, not a shared crop -- each frame's already been independently
# content-trimmed to its own extent by the render_* call that produced it,
# so a naive shared crop can't realign them). Skips a view entirely if
# fewer than 2 of its frames rendered (nothing to pad to, and padding a
# single image to its own size would just be a no-op anyway).
for view_name in top lat xray trans; do
  case $view_name in
    top) pngs=("${top_pngs[@]}") ;;
    lat) pngs=("${lat_pngs[@]}") ;;
    xray) pngs=("${xray_pngs[@]}") ;;
    trans) pngs=("${trans_pngs[@]}") ;;
  esac
  if [ "${#pngs[@]}" -lt 2 ]; then
    echo "  skipping size normalization for $view_name (${#pngs[@]} frame(s))"
    continue
  fi
  echo "Normalizing $view_name to a uniform size across ${#pngs[@]} frames"
  workspace_pngs=("${pngs[@]/#/\/workspace\/}")
  docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 -v "$(pwd)":/workspace --entrypoint /opt/paraview/bin/pvpython \
    kitware/paraview:pv-v5.8.0-osmesa-py3 \
    /workspace/results/render_view.py --trim-batch "${workspace_pngs[@]}" \
    > /tmp/stackvid_trim_${PREFIX}_${view_name}.log 2>&1 || { echo "  uniform-trim FAILED for $view_name (see /tmp/stackvid_trim_${PREFIX}_${view_name}.log)"; }
done

# Pass 3: stack. Order: top, transverse, lateral, xray. Colorbars are
# embedded inside top_png/lat_png already (render_view.py's
# _save_colorbar()), so this is just a straight 4-input vstack -- no
# separate legend row to build, and no per-panel cropping (the old
# transverse bottom-20% crop was dropped, user request, 2026-08-03: no
# longer needed now that render_transverse's own framing/margins are tuned
# directly). Every view's own frames are now a uniform size (pass 2 above),
# so W is the same on every iteration too -- kept as a per-iteration
# ffprobe computation anyway rather than hoisted out, since it's now cheap
# insurance against a view that got skipped above (single-frame batch).
i=0
for t in "${TIMES[@]}"; do
  i=$((i+1))
  top_png="results/${PREFIX}_top_screenshot_t${t}.png"
  lat_png="results/${PREFIX}_lateral_screenshot_t${t}.png"
  xray_png="results/${PREFIX}_lateral_xray_t${t}.png"
  trans_png="results/${PREFIX}_transverse_screenshot_t${t}.png"
  stacked_png="results/${PREFIX}_stacked_t${t}.png"
  if [ ! -f "$top_png" ] || [ ! -f "$lat_png" ] || [ ! -f "$xray_png" ] || [ ! -f "$trans_png" ]; then
    echo "[$i/${#TIMES[@]}] t=$t: skipping stack (a per-view render failed above)"
    continue
  fi

  echo "[$i/${#TIMES[@]}] t=$t: stacking"
  docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc "
    set -e
    W=\$(for f in /workspace/$top_png /workspace/$lat_png /workspace/$xray_png /workspace/$trans_png; do
      ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 \"\$f\"
    done | sort -n | tail -1)
    ffmpeg -y \
      -i /workspace/$top_png \
      -i /workspace/$trans_png \
      -i /workspace/$lat_png \
      -i /workspace/$xray_png \
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
# -vf scale=trunc(iw/2)*2:trunc(ih/2)*2 -- libx264 requires even width/height;
# the stacked frames' width comes from the widest per-view render (odd or
# even depending on the case's own track length), so crop the rare
# trailing odd pixel rather than let the encoder reject the whole video.
docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc \
  "ffmpeg -y -f concat -safe 0 -i /workspace/$CONCAT_LIST -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -vsync vfr -pix_fmt yuv420p -c:v libx264 -crf 18 /workspace/$VIDEO_OUT"

echo "VIDEO DONE: $VIDEO_OUT"
