#!/bin/bash
# Stacks two or more already-built cutaway_updated_video.mp4 files (from
# _render_cutaway_updated_video.sh) into one labeled, vertically-stacked
# comparison video -- used for the T0-sweep comparison (see T0_sweep.md).
#
# Each input video must already be time-synced to the same target
# timesteps (same frame count/duration) -- either because both cases use
# the same fork-stage time window/writeInterval (true for all of
# testrun75/78/81/84, so their cutaway_updated_video.mp4 line up with no
# extra work), or, for a case outside that shared grid (e.g. testrun69,
# a different lineage with irregular actual write times), by rendering
# it at the other cases' own explicit time list instead of its own
# auto-detected timesteps. Mismatched durations will visibly desync
# partway through.
#
# What this does to each panel, top to bottom:
#   - crops the bottom 20% off (trims empty build-plate substrate so the
#     stack isn't excessively tall)
#   - draws its label (top-left, white on translucent black) at fontsize
#     60 (2.5x ffmpeg's normal default text size, for legibility once
#     stacked into a tall video)
#   - vstacks all inputs top to bottom in the order given, then re-encodes
#     at a real 10fps CFR (not the source concat-demuxer VFR) with -g 10
#     (1 keyframe/sec) so seeking works correctly in any player, VLC
#     included -- see the same fix in _render_cutaway_updated_video.sh's
#     own header comment for why this matters.
#
# Usage (N >= 2 video/label pairs, top to bottom):
#   bash results/_build_cutaway_updated_collage.sh <output.mp4> \
#     <video_1> <label_1> [<video_2> <label_2> ...]
#
# Example (Zixun's 500K-over-400K comparison):
#   bash results/_build_cutaway_updated_collage.sh \
#     results/T0_comparison/collage_500_400.mp4 \
#     results/testrun84_vdep_3_Al_cutaway_updated_video.mp4 "testrun84 (T0=500K)" \
#     results/testrun81_vdep_3_Al_cutaway_updated_video.mp4 "testrun81 (T0=400K)"
#
# Example (this session's 3-way T0=100/200/300K collage):
#   bash results/_build_cutaway_updated_collage.sh \
#     results/T0_comparison/cutaway_updated_collage_tr69_tr75_tr78.mp4 \
#     results/tr69/cutaway_updated_video.mp4 "testrun69 (T0=300K)" \
#     results/tr78/cutaway_updated_video.mp4 "testrun78 (T0=200K)" \
#     results/tr75/cutaway_updated_video.mp4 "testrun75 (T0=100K)"
set -euo pipefail
cd ~/lbf3

if [ $# -lt 5 ] || [ $(( ($# - 1) % 2 )) -ne 0 ]; then
  echo "Usage: bash results/_build_cutaway_updated_collage.sh <output.mp4> <video_1> <label_1> [<video_2> <label_2> ...]"
  echo "  (at least 2 video/label pairs, in top-to-bottom order)"
  exit 1
fi

OUT="$1"; shift
VIDEOS=(); LABELS=()
while [ $# -gt 0 ]; do
  VIDEOS+=("$1"); LABELS+=("$2"); shift 2
done
N=${#VIDEOS[@]}

for f in "${VIDEOS[@]}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: input video not found: $f"
    exit 1
  fi
done

mkdir -p "$(dirname "$OUT")"

INPUT_ARGS=""
FILTER=""
VSTACK_INPUTS=""
for i in "${!VIDEOS[@]}"; do
  INPUT_ARGS="$INPUT_ARGS -i /workspace/${VIDEOS[$i]}"
  FILTER="$FILTER[$i:v]crop=iw:trunc(ih*0.8/2)*2:0:0,drawtext=text='${LABELS[$i]}':x=10:y=10:fontsize=60:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=5[v$i];"
  VSTACK_INPUTS="$VSTACK_INPUTS[v$i]"
done
FILTER="$FILTER${VSTACK_INPUTS}vstack=inputs=${N}[out]"

docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc "
ffmpeg -y $INPUT_ARGS \
  -filter_complex \"$FILTER\" \
  -map \"[out]\" -r 10 -vsync cfr -g 10 -pix_fmt yuv420p -c:v libx264 -crf 18 /workspace/$OUT
"

echo "COLLAGE DONE: $OUT"
