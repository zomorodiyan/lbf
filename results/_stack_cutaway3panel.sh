#!/bin/bash
# Stacks one timestep's top/cutaway/xray PNGs (already rendered by
# render_view.py --view={top,cutaway,xray}) into a single 3-panel
# composite: top view, "Lateral-Cutaway" (the x<-25um half-domain cut,
# see render_cutaway), "Lateral-Attenuation" (the xray view). Companion
# to render_view.py's own --view=cutaway -- everything here is
# compositing/labeling only, no rendering.
#
# Usage: _stack_cutaway3panel.sh <top.png> <cutaway.png> <xray.png> \
#          <output.png> <cutaway_colorbar.png> <time_seconds>
#
# <cutaway_colorbar.png> is the standalone colorbar render_cutaway saves
# via _overlay_colorbar(..., skip_overlay=True) -- its own per-panel
# canvas can't paste it overlapping the *adjacent* top panel, so that
# overlay happens here instead, once both panels already share one
# canvas (see the H1/OVERLAP_PX block below).
#
# All label positions/sizes/text below were tuned interactively over
# many rounds (2026-08-17) -- treat the specific numbers as the current
# answer, not self-evidently "correct" ones, if revisiting.
set -euo pipefail
TOP="$1"; CUT="$2"; XRAY="$3"; OUT="$4"; CUT_COLORBAR="$5"; TIME_S="$6"
cd /workspace
W=$(for f in "$TOP" "$CUT" "$XRAY"; do
  ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 "$f"
done | sort -n | tail -1)

# Height the TOP panel will have after being scaled to width W -- needed
# to know where the top/cutaway seam lands in the final stacked image, so
# the cutaway's standalone colorbar can be pasted straddling that seam,
# poking up into the top panel ("that moves it out of the lateral-cutaway
# view and slightly into the top view but that is cool").
TOP_W=$(ffprobe -v error -select_streams v:0 -show_entries stream=width -of csv=p=0 "$TOP")
TOP_H=$(ffprobe -v error -select_streams v:0 -show_entries stream=height -of csv=p=0 "$TOP")
H1=$(( TOP_H * W / TOP_W ))
OVERLAP_PX=40  # how far the colorbar pokes up into the top panel

# All labels top-left, black, no background box -- top-left sits against
# each panel's own white "gas" headspace region, so plain black reads
# fine everywhere, including over the xray panel's mostly-black lower area.
DT_BLACK="fontsize=40:fontcolor=black:x=70:y=40"
# Top panel's own axis-direction text, slightly larger, immediately after
# "Top view".
DT_TOPLABEL2="fontsize=48:fontcolor=black:x=260:y=40"

# Timestamp: seconds -> microseconds, 2 decimal places, fontsize=56 (was
# 80 -- reduced 30%). Bottom-left of the top panel -- fully open space
# there (the panel's own "Top view x down z right" label lives top-left).
TIME_US=$(awk "BEGIN {printf \"%.2f\", ${TIME_S} * 1e6}")
DT_TIME="fontsize=56:fontcolor=black:x=70:y=h-th-20"

ffmpeg -y -loglevel error -i "$TOP" -i "$CUT" -i "$XRAY" -i "$CUT_COLORBAR" \
  -filter_complex "
    [0:v]scale=${W}:-2,drawtext=text='Top view  ':${DT_BLACK},drawtext=text='x ↓  z →':${DT_TOPLABEL2},drawtext=text='t = ${TIME_US} us':${DT_TIME}[v0];
    [1:v]scale=${W}:-2,drawtext=text='Lateral-Cutaway  x < -25um':${DT_BLACK}[v1];
    [2:v]scale=${W}:-2,drawtext=text='Lateral-Attenuation':${DT_BLACK}[v2];
    [v0][v1][v2]vstack=inputs=3[stacked];
    [stacked][3:v]overlay=x=W-w-15:y=${H1}-${OVERLAP_PX}
  " \
  "$OUT"
echo "Stacked -> $OUT"
