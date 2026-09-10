#!/bin/bash
# Stacks three already-built cutaway_updated_video.mp4 files (from
# _render_cutaway_updated_video.sh) into one labeled, vertically-stacked
# comparison video -- used for the T0-sweep comparison (see T0_sweep.md).
#
# Each input video must already be time-synced to the same target
# timesteps (same frame count/duration) -- e.g. by rendering each case at
# the same explicit time list rather than each case's own auto-detected
# timesteps, the way testrun69's frames were rendered against tr75's own
# time list. Mismatched durations will visibly desync partway through.
#
# What this does to each panel, top to bottom:
#   - crops the bottom 20% off (trims empty build-plate substrate so the
#     stack isn't excessively tall)
#   - draws its label (top-left, white on translucent black) at fontsize
#     60 (2.5x ffmpeg's normal default text size, for legibility once
#     stacked into a tall video)
#   - vstacks all three, then re-encodes at a real 10fps CFR (not the
#     source concat-demuxer VFR) with -g 10 (1 keyframe/sec) so seeking
#     works correctly in any player, VLC included -- see the same fix in
#     _render_cutaway_updated_video.sh's own header comment for why this
#     matters.
#
# Usage:
#   bash results/_build_cutaway_updated_collage.sh <output.mp4> \
#     <video_top> <label_top> <video_mid> <label_mid> <video_bottom> <label_bottom>
#
# Example (this session's T0=100/200/300K collage):
#   bash results/_build_cutaway_updated_collage.sh \
#     results/T0_comparison/cutaway_updated_collage_tr69_tr75_tr78.mp4 \
#     results/tr69/cutaway_updated_video.mp4 "testrun69 (T0=300K)" \
#     results/tr78/cutaway_updated_video.mp4 "testrun78 (T0=200K)" \
#     results/tr75/cutaway_updated_video.mp4 "testrun75 (T0=100K)"
set -euo pipefail
cd ~/lbf3

if [ $# -ne 7 ]; then
  echo "Usage: bash results/_build_cutaway_updated_collage.sh <output.mp4> <video_top> <label_top> <video_mid> <label_mid> <video_bottom> <label_bottom>"
  exit 1
fi

OUT="$1"; V_TOP="$2"; L_TOP="$3"; V_MID="$4"; L_MID="$5"; V_BOT="$6"; L_BOT="$7"

for f in "$V_TOP" "$V_MID" "$V_BOT"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: input video not found: $f"
    exit 1
  fi
done

mkdir -p "$(dirname "$OUT")"

docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)":/workspace lbf3 bash -lc "
ffmpeg -y \
  -i /workspace/$V_TOP \
  -i /workspace/$V_MID \
  -i /workspace/$V_BOT \
  -filter_complex \"
    [0:v]crop=iw:trunc(ih*0.8/2)*2:0:0,drawtext=text='$L_TOP':x=10:y=10:fontsize=60:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=5[v0];
    [1:v]crop=iw:trunc(ih*0.8/2)*2:0:0,drawtext=text='$L_MID':x=10:y=10:fontsize=60:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=5[v1];
    [2:v]crop=iw:trunc(ih*0.8/2)*2:0:0,drawtext=text='$L_BOT':x=10:y=10:fontsize=60:fontcolor=white:box=1:boxcolor=black@0.5:boxborderw=5[v2];
    [v0][v1][v2]vstack=inputs=3[out]
  \" \
  -map \"[out]\" -r 10 -vsync cfr -g 10 -pix_fmt yuv420p -c:v libx264 -crf 18 /workspace/$OUT
"

echo "COLLAGE DONE: $OUT"
