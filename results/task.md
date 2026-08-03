# results/ post-processing task list

Working list of improvements to the post-processing pipeline --
`render_view.py` (the four views: top, lateral, xray, transverse, formerly
4 separate scripts, merged -- see Done below) and `_render_stacked_video.sh`.
Add tasks here as they come up; move to Done as they're addressed.

## Open

### 1. Trim empty borders / tighten crop ranges to shrink image sizes

Cut off unnecessary blank space so the output images (and the stacked
composite) are smaller in both pixel dimensions and file size.

**Vertical whitespace: solved** (see Done below) -- `_trim_vertical_whitespace()`
(top/lateral) and `_trim_vertical_whitespace_uniform()` (transverse's 3
panels together) now crop the blank top/bottom margin `FRAME_MARGIN`
leaves. Unlike the horizontal case below, this needed no min-run-length
trick: the margin sits entirely *outside* the crop box, genuinely empty,
not a thin sliver of real content.

**Horizontal: still not fully solved, and now better understood.**
Investigated and ruled out two of the originally-planned sub-approaches
(both wrong once actually checked against the render pipeline):

- **Shrinking `FRAME_MARGIN` (1.3 -&gt; 1.05) backfires.** The plan was: now
  that colorbar separation (Done) removed that margin's original
  justification, shrink it to tighten the frame. Tried on
  `lateral_screenshot.py` and `top_screenshot.py`, measured on testrun64 at
  t=1.974e-4s: post-trim kept-column width went from 2503px -&gt; 3094px
  (lateral) and 2443px -&gt; 3017px (top) -- both *larger* files, not smaller.
  Root cause: zooming in (smaller margin) makes the untouched flat plate's
  edge-on silhouette taller in *pixels*, and the existing post-render
  side-trim tells that silhouette apart from real content via a *fixed*
  `_MIN_COL_HEIGHT_PX = 5` threshold -- so a smaller margin pushes more of
  the "boring" background past that fixed threshold, and the trim ends up
  keeping *more* columns. Reverted to 1.3 in both scripts (with the finding
  documented inline).
- **`transverse_screenshot.py` doesn't need a post-render trim step at
  all.** It was flagged as missing one, but its own header already explains
  why it doesn't need one: unlike `lateral_screenshot.py`/`top_screenshot.py`
  (whose z-window spans the *entire* scan, so short/early tracks leave
  genuine blank space to trim), `transverse_screenshot.py`'s per-panel crop
  window is small and fixed around each cut, not scan-spanning -- there's no
  accumulated blank space in the first place. (Also fixed a stale comment
  found along the way: `FRAME_MARGIN`'s justification there referenced "the
  colorbar/text overlay," but this script has never had either in the
  current codebase -- see its own "No colorbar/text label overlay" comment.
  Left the value at 1.3 regardless, same reasoning as above.)
- **`lateral_xray.py` likely doesn't need one either, for a different
  reason**: its output isn't a 3D render with white background -- every
  pixel is real Beer-Lambert attenuation data (matplotlib grayscale), so the
  "wasted space" for an early/short track isn't blank at all, it's a
  genuinely-correct dark region representing solid, unmelted plate. There's
  no empty margin there to trim; the fixed-z-window "problem" the other two
  scripts have doesn't apply to this script's output in the same way. Not
  yet independently re-verified by rendering (unlike the two points above,
  which were), so treat this one as a hypothesis pending a quick check
  before ruling it out for good.

**New finding while adding vertical trim**: `top_screenshot`'s (now
`render_top`'s) horizontal trim doesn't actually remove much blank space
for early/short tracks, and this is structural, not a bug. The
`_MIN_COL_HEIGHT_PX=5` min-run trick exists to ignore the untouched flat
plate's *edge-on, ~1px-tall* silhouette in `render_lateral` -- but
`render_top` looks straight down *at* the plate, not edge-on, so its
untouched-plate background is a genuinely tall (not thin-sliver) band
across the whole crop width, easily clearing that 5px threshold and
counting as "real content" everywhere. This isn't wasted canvas exactly
(it's a real, correctly-rendered surface, just visually plain/undisturbed)
but it does mean `render_top` frames don't shrink for early timesteps the
way `render_lateral`'s do. Confirmed pre-existing (present before today's
changes too), not a regression.

**What's actually left to investigate**: tightening the *physical*
`Y_DEPTH_MIN/MAX`/`X_LATERAL_MIN/MAX` window constants themselves (not
camera zoom) -- these are fixed, hand-picked values, not measured from the
data, and directly set the rendered canvas dimensions (via the aspect-ratio
computation feeding `ViewSize`), so shrinking them for real would shrink the
canvas itself rather than interacting with the trim heuristic the way
`FRAME_MARGIN` did. This needs checking across cases/timesteps (not just
testrun64) to make sure the window isn't cut tighter than the melt-affected
region ever needs -- higher risk of clipping real content than the two
ruled-out approaches, so worth confirming direction before implementing.

Notes:
- The stacked composite (`_render_stacked_video.sh`) vstacks the 4 views
  (colorbars now embedded inside top/lateral's own images, no separate
  legend row -- see Done) -- any leftover per-view margin compounds
  there, so tightening the individual views should shrink the
  composite/mp4 too.

## Done

### Schematic follow-ups, transverse depth clip + colored cap sides, narrower X_LATERAL

Five more changes, from the round of feedback right after the previous
entry below:

1. **Schematic "z" label nudged 4 characters further down**, close to the
   bottom of the domain box -- a screen-space pixel offset via
   `ScaledTranslation` composed onto the label's own transform (same
   technique the removed tick-label nudging used, minus the "reapply every
   draw" patch that was only needed for mplot3d-managed tick Text objects,
   not a plain `ax3d.text()` artist).
2. **Schematic rotated back slightly**: `azim` `-20 -> -30` (partway back
   toward the original `-55`, user request "slightly retotate back" after
   the previous round's `-20` turned out a bit too flat).
3. **Transverse panels clipped to 200um depth** (was the shared 400um) --
   local `Y_DEPTH_MAX = 0.2e-3` override inside `render_transverse` only
   (shadows the module-level constant used by top/lateral, same pattern as
   `render_xray`'s own `XRAY_Y_DEPTH_MAX`), affecting the `box_clip`
   window, `y_center`, and `view_w`.
4. **Cross-section cap's straight sides colored to match each panel's own
   rim color** (previously plain gray, same as the cap's fill) -- 3
   explicit `Line` sources per panel (left/right/bottom of the box-clip
   crop window) in that panel's `gm_rim_color`, rather than enabling
   `EdgeVisibility` on the cap surface itself (which would draw every
   individual mesh cell's edges, not just this outer silhouette).
   `gm_rim_color` computation moved earlier in the per-panel loop so both
   this and the existing `gm_rim_disp` can share it.
5. **`X_LATERAL_MIN`/`MAX` narrowed from +/-0.2mm to +/-0.18mm** -- shared
   between top and transverse deliberately (confirmed with the user first
   that narrowing both together, not just top, was intended).
6. **Transverse's bottom 10% cropped before stacking**, applied only in
   `_render_stacked_video.sh`'s ffmpeg filter (`crop=iw:ih*0.9:0:0` on the
   transverse input before its `scale`), not to the standalone
   `demo_*_transverse.png`/production per-frame transverse file itself.

Verified each piece individually (transverse showing visibly less depth
with colored cap sides; top view's narrower vertical extent; schematic's
z-label position and rotation) and the full 4-timestep demo set together,
including the new stacked-composite crop.

### Transverse Y-depth window shift + cap fill/outline color swap

Two more changes to `render_transverse` on top of the previous entry:

1. **Y-depth window shifted from 50-200um to 100-500um below the surface**
   -- previously only `Y_DEPTH_MAX` was a transverse-local override
   (shadowing the shared module-level constant used by top/lateral);
   `Y_DEPTH_MIN` is now locally overridden too (`0.1e-3`), alongside
   `Y_DEPTH_MAX` (`0.5e-3`), so the window no longer starts near the
   surface and now reaches deep enough to show the full melt-pool depth
   including the keyhole bottom.
2. **Cross-section cap's fill and outline colors swapped**: the filled
   solid portion (`cap_solid_disp`) now uses the panel's own
   `gm_rim_color` (red/orange/yellow) instead of the old flat gray
   `CAP_SOLID_COLOR` (removed, now unused); the crop-boundary side
   outlines (the 3 `Line` sources added in the previous entry) now use
   plain white instead of `gm_rim_color`, so the outline stays visually
   distinct now that the fill itself carries the panel color.

Verified via a test render (`testrun64`, t=1.734e-4s) showing the deeper
window (full melt-pool bottom visible) and the swapped colors (colored
fill, white outline), then confirmed at t=3.999e-4s in the regenerated
full 4-timestep demo set.

### Schematic: gray solid surface, slight further rotation

`plot_domain_schematic.py`: metal slab `facecolors` `"silver" -> "gray"`
(user request "turn the solid surface color to gray"); `azim` `-30 -> -35`,
continuing the same "rotate back toward -55" direction as the previous
round's `-20 -> -30` step (user asked for a further slight rotation without
specifying direction -- picked the same direction as last time's step
rather than reversing it; easy to flip if wrong). Verified via a standalone
test render (no simulation data needed, script has CLI defaults).

### Transverse: cap-solid fill reverted to gray

`render_transverse`'s cross-section cap fill (`cap_solid_disp`) reverted
from the panel's own `gm_rim_color` back to flat gray (`CAP_SOLID_COLOR`,
re-added) -- user request "color the solid surface as gray like before",
matching the schematic's own gray solid slab. The crop-boundary side
outlines stay white (unchanged from the previous round); the gm_rim
outline at the actual gas/metal cut still uses each panel's red/orange/
yellow color, so the panels remain distinguishable by that outline alone.
Verified via a test render, then regenerated the full 4-timestep demo set.

### Cross-section marker palette swap + transverse Y-depth to 100-450um

Two more changes:

1. **`GM_RIM_COLORS` (the shared red/orange/yellow palette used by the
   vertical cut-position lines in top/lateral and the gm_rim outlines in
   transverse) changed to cyan/pistachio-green/yellow** -- red `[1,0,0] ->`
   cyan `[0,1,1]`, orange `[1,0.5,0] ->` pistachio green `[0.576,0.773,
   0.447]`, yellow unchanged. Single shared constant, so this one edit
   propagates to all three views without touching per-view code.
2. **Transverse Y-depth window narrowed from 100-500um to 100-450um**
   (`Y_DEPTH_MAX` local override `0.5e-3 -> 0.45e-3`).

Verified via test renders of top and transverse (all 3 colors visible and
correctly ordered, tighter depth window), then regenerated the full
4-timestep demo set (top/lateral/transverse; xray untouched since it
doesn't use `GM_RIM_COLORS`).

### Pistachio green swapped for a more saturated true green

User feedback right after the previous entry: "the green does not work but
yellow and cyan are alright" -- pistachio green `[0.576,0.773,0.447]`
(desaturated, close in luminance to the render's own orange/tan tones)
didn't read clearly against the data. Replaced with a saturated true green
`[0.0,0.8,0.0]`. Verified via a test render (clearly distinct from both
cyan and yellow, and from the render's own orange/purple/gray coloring),
then regenerated the full 4-timestep demo set.

### Transverse colorbar clipped at the bottom, lateral top-10% clip

Two more fixes:

1. **Transverse's `z - z_laser` colorbar was clipped at the bottom of the
   image** -- `_overlay_colorbar`'s default `bottom_margin=15` assumes some
   blank canvas below the content (true for top/lateral), but transverse's
   3 panels are trimmed tight with no blank margin at the bottom (the crop
   window's own flat bottom edge sits right at the image edge), so the bar
   ended up sitting flush against it. Raised `bottom_margin` to `60` for
   `render_transverse`'s `_overlay_colorbar` call only (top/lateral
   untouched, still using the default).
2. **Lateral view: clip the top 10%** -- new `_clip_top_fraction(output_png,
   frac)` helper (crops off the top `frac` of the already-trimmed image),
   called in `render_lateral` right after `_trim_vertical_whitespace`, cuts
   off excess gas headspace above the melt track/beads without touching any
   real content.

Verified both via test renders, then regenerated the full 4-timestep demo
set (lateral + transverse; top/xray untouched).

### Schematic: another slight rotation, per-axis label nudges

1. **`azim` `-35 -> -40`**, another slight rotation in the same direction
   as the last two rounds' steps.
2. **x/y/z label positions nudged**, in screen-space character units (same
   `ScaledTranslation` technique as the earlier z-nudge, deliberately
   screen-space so the offset stays fixed on the page regardless of any
   later `azim` change -- user specified these "as if there was not
   rotation"): y 1 character left, x 1 character down, z 5 characters
   right (in addition to the existing 4-characters-down nudge from a
   previous round). Refactored the per-axis nudge from a single
   z-only `if` block into a shared `{"x":..., "y":..., "z":...}` dict
   inside the label loop.

Verified via a standalone test render (no simulation data needed).

### Schematic: thicker lines, black track extension removed, further z/time-label adjustments

Follow-up round, mostly mid-turn adjustments after seeing the render:

1. **Thicker lines**: domain wireframe `1.6 -> 2.5`, red active-track line
   `3 -> 4.5`, the (now-removed, see #2) black extension lines were `2 ->
   3` before being dropped entirely.
2. **Removed the black laser-track extension lines** beyond the active red
   portion, out to the domain's z edges -- user request, "no need to
   extend the laser track all the way."
3. **Time label moved further from the beam**: its y-offset above the
   beam's top `35 -> 70um`.
4. **z label pushed further right**: `5 -> 9 -> 12` characters (two
   follow-up nudges after the initial `+5`), same screen-space nudge
   pattern as the other axis labels.

Verified via standalone test renders after each change, then regenerated
the full 4-timestep demo set (transverse only, since the schematic is
embedded per-frame there; top/lateral/xray untouched).

### Schematic axis labels: bigger + direction arrows; colorbars moved bottom-right

Two more changes:

1. **Schematic x/y/z labels bigger and carry a direction arrow**: new
   `AXIS_FONTSIZE = round(TIME_FONTSIZE * 1.3)` (was sharing `TIME_FONTSIZE`
   with the time label) and each label's text changed from the bare letter
   to `"x↘"`, `"y↓"`, `"z↗"` (plain Unicode arrow glyphs -- DejaVu Sans,
   this Docker image's default font, covers them fine, unlike the missing
   micro-sign glyph noted elsewhere in this repo). No space between the
   letter and arrow (user request). `y`'s arrow was originally `↑`,
   corrected to `↓` per user follow-up.
2. **All 3 colorbars (top/lateral/transverse) moved from bottom-center to
   bottom-right** of their respective images -- `_overlay_colorbar` gained
   a `right_margin=15` parameter; the overlay's `x0` is now computed as
   `vw - right_margin - cb_img.shape[1]` instead of centering. Shared
   helper, so this one change applies to all three views at once.

Verified via test renders of the schematic and top/transverse, then
regenerated the full 4-timestep demo set (top/lateral/transverse; xray
untouched, it never had a colorbar).

### Colorbar sizing made proportional (cross-view consistency), schematic fonts bumped further, time label styling

User feedback: colorbar/text sizes weren't consistent across the stacked
rows, and wanted the schematic's axis-label fonts even bigger plus a
different look for the time label.

1. **`_overlay_colorbar`'s `cb_height` (and its font sizes) now scale with
   `output_png`'s own width**, instead of a fixed `cb_height=100`.
   `_render_stacked_video.sh` later scales every view's saved PNG to one
   shared width for the stacked composite -- since top/lateral/transverse
   each have a different native width, a fixed absolute `cb_height`
   ended up a different *effective* size post-scale for each. Anchoring it
   to `vw` (`cb_height = round(vw * 100/1633)`, `1633` = top view's own
   historical width, kept as the reference so top's own look is
   unchanged) makes the eventual post-scale size depend only on the
   shared target width, not each view's native resolution.
2. **`render_transverse`'s colorbar overlay moved to *after* the schematic
   prepend** (was before) -- the proportional sizing above only works if
   the width it measures is the row's true final width; previously the
   colorbar was sized against the pre-schematic panel-only composite
   (~1260px) while the actually-saved file was wider (~1839px,
   schematic included), undermining the fix above for transverse
   specifically.
3. **Schematic axis-label font bumped further**: `AXIS_FONTSIZE`
   multiplier `1.3 -> 1.6` of `TIME_FONTSIZE`.
4. **Time label restyled**: color `darkorange -> black`, size
   `TIME_FONTSIZE -> round(TIME_FONTSIZE * 1.15)`, offset above the beam
   raised again `70 -> 100um` ("a bit higher").

Verified the colorbar fix by simulating the stacking script's own
scale-to-common-width step on freshly rendered top/lateral/transverse
files and comparing cropped bottom-right corners side by side (text sizes
now match); verified the schematic changes via a standalone test render.
Regenerated the full 4-timestep demo set (top/lateral/transverse; xray
untouched).

### Colorbars unified to mm, tighter transverse gaps

1. **Top/lateral colorbars converted from um to mm**, matching
   transverse's `z - z_laser` (already mm) -- `render_top`'s `ycolor`
   Calculator (`*1e6 -> *1e3`), `ymin_off`/`ymax_off`, `transition`
   (100um -> 0.1mm), title `'y (um)' -> 'y (mm)'`, `custom_labels`
   `[-100,0,100] -> [-0.1,0,0.1]`; same changes mirrored in `render_lateral`
   for `xcolor`/`xmin_mm`/`xmax_mm`/`'x (mm)'`.
2. **Transverse colorbar title changed** `'z - z_laser' -> 'z - z_l (mm)'`
   (its values were already mm, just the title lacked a unit).
3. **Tighter transverse layout**: `PANEL_GAP_PX` `6 -> 3`,
   `SCHEMATIC_GAP_PX` `20 -> 8`, and `_trim_side_whitespace_uniform`'s own
   per-panel side padding `10 -> 3` (this padding applies on *both* sides
   of *every* panel, so it was compounding with `PANEL_GAP_PX` between
   adjacent panels -- the dominant source of the visible gap, not
   `PANEL_GAP_PX` alone).

Verified via test renders (top: mm labels correct; transverse: panels
visibly tighter together, schematic closer to the first panel, colorbar
title updated), then regenerated the full 4-timestep demo set
(top/lateral/transverse; xray untouched, no colorbar).

### Schematic cross-section cut lines, top-view laser marker, schematic cleanup

1. **Schematic cross-section cut lines**: `plot_domain_schematic.py` now
   duplicates `OFFSETS_BEHIND_LASER`/`GM_RIM_COLORS` (same "small constants
   duplicated per standalone pvpython script" convention as
   `top_screenshot.py`) and, for each cut whose `z = LASER_Z_UM - offset`
   falls inside the domain, draws an L-shaped line in that cut's color:
   along the plate's top surface (spanning the full x width) and down the
   visible `x=X1` (0.32mm) side face (spanning the plate's solid depth) --
   user request, 2026-08-02 ("add cross-section lines to the schematic
   with correct distance behind z_laser on the plate's top surface and the
   x=0.32 surface").
2. **Top-view laser center-position marker**: `render_top` now shows a
   small solid orange `Sphere` (20um radius) at the laser's actual current
   `(x=0, z=laser_z)` -- x=0 because every view in this pipeline already
   assumes the laser travels along the x=0 centerline (no separate x
   lookup exists). Same orange as the schematic's own beam/landing marker
   for visual consistency. User request, 2026-08-02 ("show a simple
   appropriate sign at the laser center position in the top view").
3. **Schematic cleanup** (mid-turn follow-up after seeing the render):
   removed the translucent gas-headspace fill above the surface (now plain
   white background inside the wireframe box) and doubled the domain
   wireframe's line width again (`2.5 -> 5.0`).

Verified via standalone schematic test renders and a top-view test render
(t=1.734e-4s, where the laser sits inside the fixed Z_VIEW crop window --
t=3.999e-4s's laser position falls outside that window, so the marker
correctly doesn't appear there, same expected behavior as the existing
offset-marker lines). Regenerated the full 4-timestep demo set (top +
transverse; lateral/xray untouched).

### Schematic plate-top-rim edge restored, schematic enlarged in the transverse composite

1. **Fixed the missing plate-top-rim edge**: removing the gas-headspace
   fill (previous round) took away the only thing that had ever implied
   that edge (color contrast against the translucent fill) -- there was no
   actual line there. Fixed with `draw_box_wireframe(X0, X1, Y_SURF,
   Y_SURF, Z0, Z1, ...)`: passing `y0=y1=Y_SURF` collapses
   `box_vertices()` into a degenerate zero-height box, so the 4 "vertical"
   edges vanish (zero-length) and exactly the 4 perimeter edges of the
   y=Y_SURF rectangle remain, each still going through the existing
   occlusion logic (near/right edges visible, far/left hidden) -- reusing
   the tested function instead of writing a one-off 4-edge helper.
2. **Schematic enlarged within the transverse composite**: was resized to
   exactly match the panel row's own height (1:1), reading too small/hard
   to make out (user report, 2026-08-02). Now resized to
   `row_img.shape[0] * 1.35` (`SCHEMATIC_HEIGHT_FACTOR`), with the panel
   row white-padded top/bottom (centered) to match the taller schematic
   instead of the other way around. `SCHEMATIC_GAP_PX` narrowed further
   (8 -> 4) to partly offset the extra width from the larger schematic, per
   the user's own suggestion.

Verified via a standalone schematic test render (plate edge now visible on
both the front and right/x=X1 faces) and a transverse test render (visibly
larger, more legible schematic with white-padded panels alongside it),
then regenerated the full 4-timestep demo set (transverse only;
top/lateral/xray untouched).

### Transverse: colorbar compensated for the stacked-composite crop, doubled panel gap, no schematic gap

1. **Colorbar bottom margin now a fraction of image height, not a fixed
   pixel count**: `_overlay_colorbar` gained `bottom_margin_frac` (used
   instead of `bottom_margin` when given); `render_transverse` now passes
   `bottom_margin_frac=0.18`. `_render_stacked_video.sh` crops the bottom
   10% off the transverse row for the *stacked composite only* (not this
   standalone file) -- a fixed pixel margin cleared that crop line by a
   different, sometimes too-small amount depending on each frame's own
   image height; a fraction-based margin clears it by the same proportion
   every time. Verified by simulating that same 10% crop on a fresh
   render and confirming the bar stays fully intact with margin to spare.
2. **Panel-to-panel gap doubled**: `PANEL_GAP_PX` `3 -> 6` and
   `_trim_side_whitespace_uniform`'s own per-panel side pad `3 -> 6`
   (scaled together since the total visible gap is `2*pad + PANEL_GAP_PX`)
   -- user request, 2026-08-02, "the space in between transverse images
   should [be] twice what it is."
3. **No gap between the schematic and the first panel**: `SCHEMATIC_GAP_PX`
   `4 -> 0`.

Regenerated the full 4-timestep demo set (transverse only; top/lateral/
xray untouched).

### More transverse crop + colorbar compensation, bigger time font, smaller xray bar text, dark-gray plate rim

1. **`_render_stacked_video.sh`'s bottom crop on the transverse row
   doubled**: `ih*0.9` (10% cut) `-> ih*0.8` (20% cut) -- user request,
   2026-08-02, "we need more cut-off from bottom." Applied only to the
   stacked composite, not the standalone `demo_*_transverse.png`.
   `render_transverse`'s colorbar `bottom_margin_frac` raised to match:
   `0.18 -> 0.28`, keeping the same ~0.08 buffer past the crop line as
   before. Verified by simulating the new 20% crop on a fresh render --
   the bar stays fully intact with margin to spare.
2. **Schematic time-label font multiplied by another 1.5x**: net
   multiplier on `TIME_FONTSIZE` is now `1.15 * 1.5 = 1.725` (was `1.15`).
3. **Xray scale bar's "0.5mm" text shrunk to 0.75x**: `BAR_FONTSIZE`
   `20 -> round(20 * 0.75) = 15`. Clarified with the user first that
   "remove" in their message meant resize, not delete (the rest of the
   message was about font-size scaling, and "remove ... to 0.75x" is
   self-contradictory as a literal deletion).
4. **Schematic plate-top-rim recolored dark gray** (was black, same as the
   true domain-boundary wireframe): `draw_box_wireframe(X0, X1, Y_SURF,
   Y_SURF, Z0, Z1, color="dimgray", ...)` -- user request, mid-turn
   follow-up after seeing the previous render, since this rim is a plate
   edge, not an actual domain-boundary edge, and should read as visually
   distinct from the edges that are.

Verified each individually via test renders, then regenerated the full
4-timestep demo set (xray + transverse; top/lateral untouched).

### Top-view laser marker: hollow circle at the true spot radius, not a filled sphere

User request: "draw a circle with the right radius for the laser and a
dot ... in the middle ... in orange [with white for the dot]" -- replaced
`render_top`'s plain filled orange `Sphere` marker with two pieces: a
hollow orange circle at the laser's true physical spot radius (35um, see
CLAUDE.md's "Physical parameters" section) plus a small white dot at its
center. No `RegularPolygonSource`/`Circle`/annulus-outline source exists
in this ParaView build (checked `dir(paraview.simple)` directly -- only
`Disk`, which would give a filled/banded shape, not a clean outline), so
the circle is built by hand via `PolyLineSource`: 48 points around the
circle in the x/z plane (`Closed=1` connects the last point back to the
first), `Show(..., Representation='Wireframe')`. The center dot stays a
small `Sphere` (radius 6um, white, was the marker itself at 20um orange).

Verified via a test render (t=1.734e-4s), then regenerated the full
4-timestep demo set (top only; transverse/lateral/xray untouched).

### Schematic plate-rim darkened, transverse cap lightened, colorbar tick labels moved below the bar

1. **Schematic plate-top-rim darkened**: `"dimgray" -> (0.25,0.25,0.25)`.
2. **Transverse cross-section cap fill lightened**: `CAP_SOLID_COLOR`
   `[0.8,0.8,0.8] -> [0.9,0.9,0.9]`.
3. **Colorbar tick value labels moved below the bar** (were above, the
   ParaView default): `colorbar.TextPosition = 'Ticks left/bottom,
   annotations right/top'` in `_overlay_colorbar` -- found via
   `cb.GetProperty('TextPosition').GetDomain('enum')` (no such option was
   otherwise documented/obvious; the two enum values are "Ticks
   right/top, annotations left/bottom" (default) and this one). Applies to
   all three colorbar-bearing views (top/lateral/transverse) via the
   shared helper. `_find_bar_row_range`'s docstring updated to match (the
   function itself was already order-independent, just mislabeled).

Verified via test renders (schematic: visibly darker rim; transverse:
lighter cap fill and labels now below the bar), then regenerated the full
4-timestep demo set (top/lateral/transverse; xray untouched, no
colorbar).

### Schematic track line split red/dotted-white at the laser, top-view marker's center dot removed

1. **Schematic laser track line split at the laser's current position**:
   red from `TRACK_Z0` up to the laser (was red for the whole
   `TRACK_Z0`-`TRACK_Z1` active region regardless of where the laser
   actually is), dotted white for the remaining portion ahead of it up to
   `TRACK_Z1` -- user request, 2026-08-02. Clamped to
   `[TRACK_Z0,TRACK_Z1]` for the split point in case the laser hasn't
   started yet or has already passed the end of the active region.
2. **Top-view laser marker's white center dot removed** -- the orange
   circle (added last round) reads fine on its own; the dot was extra
   (user request, 2026-08-02: "the orange circle suffice[s]").

Verified via test renders (schematic: track visibly red-then-dotted-white
at the laser; top: circle only, no dot), then regenerated the full
4-timestep demo set (top + transverse; lateral/xray untouched).

### Transverse colorbar nudged down 2.5 characters

`_overlay_colorbar` gained a `down_shift_chars` parameter -- same "0.6 *
fontsize" screen-space character-height convention used elsewhere in this
repo (e.g. plot_domain_schematic.py's axis-label nudges), computed from the
colorbar's own (already font_scale-adjusted) `LabelFontSize`, captured
into a local variable before `Delete(cb_view)` for clarity. Positive values
reduce the block's own clearance from the bottom edge, i.e. move it down.
`render_transverse` now passes `down_shift_chars=2.5` (user request,
2026-08-02) alongside its existing `bottom_margin_frac=0.28`.

Verified via a test render, then re-simulated `_render_stacked_video.sh`'s
own 20% bottom crop on it to confirm the colorbar (and its -1.7/-1.2/-0.7
labels) still stays fully intact, just closer to the edge than before.
Regenerated the full 4-timestep demo set (transverse only; top/lateral/
xray untouched).

### Dynamic per-frame schematic, marker-excluded top/lateral crop, xray scale bar nudge

Several related changes in one round:

1. **xray's scale bar nudged down** by half a character's height (computed
   from the label's own fontsize and the figure's actual px/mm scale, not a
   hardcoded offset, so it stays correct if either changes).
2. **`plot_domain_schematic.py` is now regenerated per-frame** instead of
   being a static pre-baked image -- its single laser marker/time label
   must match the actual timestep being rendered. `render_transverse` now
   invokes it as a subprocess (`/opt/paraview/bin/pvpython
   plot_domain_schematic.py <laser_z_um> <time_label> <output_png>` --
   `sys.executable` resolves to an internal VTK launcher inside pvpython,
   not a usable interpreter, so the same hardcoded pvpython path every
   shell script in this repo already uses is used here too) with this
   frame's own `laser_z` (already computed earlier in `render_transverse`
   for the offset markers) and a `"<N> us"` time label. Skips the schematic
   (logs and continues) if the subprocess fails, rather than failing the
   whole batch render.
3. **`plot_domain_schematic.py` itself heavily simplified** (user request):
   removed the two refinement-box wireframes and the legend entirely;
   removed all numeric tick labels on every axis (`set_xticks([])` etc.)
   and matplotlib's own axis spine lines (`axis.line.set_visible(False)` --
   these are drawn independently of tick labels and, once our custom
   padded axis limits were introduced for the box-aspect workaround below,
   would otherwise visibly jut out past the hand-drawn domain wireframe);
   `x`/`y`/`z` axis labels are now placed by hand (`ax3d.text`, near each
   axis's own edge midpoint pushed a bit further outward from the box
   center) rather than via `set_xlabel`/`set_ylabel`/`set_zlabel`, both
   because `fontsize=` passed to those is silently ignored by this Docker
   image's matplotlib 3.1.1 (`.set_fontsize()` on the label object directly
   works reliably instead) and because the user wanted them "very close to
   the domain boundaries," which matplotlib's own auto-placement didn't
   reliably give even with reduced `labelpad`; changed from 3 fixed
   illustrative timestamps (0/100/400 us) to the single actual current
   laser position, with its landing marker switched from a small flat
   disk (now-unused `LASER_R`) to a fixed-pixel-size dot (so its size
   doesn't depend on the view's data scale) and its beam line/marker made
   fully solid and thicker (`alpha=1.0`, `linewidth=6->10`); the fixed
   view angle rotated from `azim=-55` to `azim=-20` ("more left to right"
   -- makes the laser-track line render nearly horizontal instead of
   diagonal; `-70`/`-90` were tried first and went the wrong direction,
   flattening the box end-on instead). A brief attempt at adding direction
   arrows to the axes was tried and dropped as not worth the trouble in
   this old matplotlib version (see git history if wanted later).
4. **Box-aspect workaround for matplotlib 3.1.1** (this Docker image's
   pvpython doesn't have `Axes3D.set_box_aspect`, added in 3.3): confirmed
   by testing that without it, mplot3d maps *every* axis's own view-limit
   span to the same on-screen cube edge regardless of true data range,
   rendering as a plain cube instead of the intended elongated-in-z
   proportions (user-flagged: "don't change the dimensions of the cube").
   Reproduced the same visual proportions by *padding* the two shorter
   axes' view limits (`set_xlim3d`/`set_ylim3d`/`set_zlim3d`, not touching
   any plotted data) so the true data occupies a smaller fraction of their
   own equal-sized cube edge, in exactly the ratio `BOX_ASPECT` specifies --
   general closed-form derivation in the code's own comment, not
   hardcoded to this domain's specific numbers.
5. **Top/lateral's horizontal crop range no longer affected by the
   transverse-cut marker lines** (user-reported: "affecting the range that
   is viewed ... not desirable"). Those lines span the full crop height by
   construction, so they trivially satisfied `_trim_side_whitespace`'s
   min-column-content-height check regardless of where they landed, pulling
   the kept-column range out to wherever a marker sat even when it was well
   past the real track's own extent. `_draw_offset_markers` now returns its
   marker `Show()` displays; new `_trim_side_whitespace_excluding_markers`
   hides them, renders once to measure content bounds, shows them again,
   renders the real (marker-inclusive) output, and crops it using the
   marker-free bounds. New `_content_col_bounds` factors the bounds-only
   half of `_trim_side_whitespace`'s logic out so both functions share it.

Verified each piece individually (schematic standalone at various rotations
before settling on the fix; render_transverse end-to-end including the
subprocess call; top at an early timestep with markers now correctly
excluded from the crop) and the full 4-timestep demo set together at the
end.

### Actually close the transverse panel gaps; white in-image scale bar for xray

Two follow-ups from the schematic-prepend round just below, once the result
was visible in practice:

1. **The 3 transverse panels still weren't touching even at
   `PANEL_GAP_PX=6`.** Root cause: each panel's own camera framing
   (`FRAME_MARGIN`) leaves blank margin on *both* left and right sides, same
   mechanism already known from top/lateral's horizontal trim -- but
   `render_transverse` had never trimmed that axis, only vertical (the
   panels being "already the same size" was true but irrelevant to the
   visible gap). Added `_trim_side_whitespace_uniform()` (mirrors
   `_trim_vertical_whitespace_uniform`: union of content columns across all
   3 panels, so they stay equal width) and call it right after the existing
   vertical trim. Panels now sit genuinely tight against the small gap.
2. **xray view: removed the numeric z-axis tick labels entirely**
   (`ax.tick_params(axis='x', bottom=False, labelbottom=False)`, ticks and
   numbers both gone now, not just the title removed in an earlier round)
   and replaced them with a white scale bar drawn *on* the image itself,
   overlaying the solid-black bottom section (`zorder=5`, safely below the
   faint 0.4mm depth-reference line) -- a bracket-style bar (`|——|`) spanning
   z=[1.5, 2.0]mm with a "0.5mm" label to its left. (User initially described
   the bar's span as 1.5-2.0mm but labeled it "1mm" -- clarified via a
   quick question that the span is correct and the label should match its
   true 0.5mm length instead.)

Verified both on 2-3 sample timesteps and the full stacked composite;
regenerated the full 4-timestep demo set.

### Prepend the domain schematic to the transverse row, tighten the panel gap

User moved `plot_domain_schematic.py`/its output `domain_schematic.png` (a
standalone, hand-built 3D reference diagram of the shared domain/refinement
boxes/laser track -- not rendered per-timestep) into `results/` and asked
for a modified version of it to sit to the left of `render_transverse`'s 3
cut panels in the composite, with less space between the panels themselves.
Clarified via a quick question: modification = crop the schematic's excess
white margins only, keep the legend and everything else as-is.

- `PANEL_GAP_PX` (the gap between the 3 cut panels) narrowed from 20 to 6px.
  New `SCHEMATIC_GAP_PX = 20` gives the schematic its own, wider separation
  from the panels so it still reads as a distinct element rather than a 4th
  cut.
- New `_trim_whitespace_bbox()`: crops an image to the bounding box of its
  non-white content on all 4 sides at once (the schematic's source figure
  has wide margins on every side) -- unlike the existing
  `_trim_side_whitespace`/`_trim_vertical_whitespace(_uniform)` helpers,
  which each only trim one axis for the ParaView-rendered views.
- New `_resize_image()`: scales an image to an exact target height,
  preserving its own aspect ratio, via the same "render through a
  throwaway matplotlib figure, read back the pixel buffer" technique
  `_render_title_image()` already uses (no PIL/scipy available in this
  Docker image for a direct resize call).
- In `render_transverse()`: after the existing colorbar overlay (so the
  colorbar stays centered under the 3 panels, not the whole row including
  the schematic), load `domain_schematic.png` (path resolved relative to
  the script's own directory, not CWD), crop and resize it to the row's
  own height, and concatenate it to the left with `SCHEMATIC_GAP_PX`
  separation. Skips cleanly (logs and continues) if the schematic file
  isn't present, e.g. if this script is later run somewhere it wasn't
  copied to.

Verified on 2 sample timesteps (t1, t2): schematic renders legibly at the
row's reduced height, panels visibly tighter, colorbar still correctly
centered under the panels. Regenerated the full 4-timestep demo set --
since transverse became the widest of the 4 stacked rows after adding the
schematic, `_render_stacked_video.sh`'s existing "scale to widest row"
logic scaled the other 3 rows up to match, unchanged code path, confirmed
still working correctly.

### Narrow the shared z-window from 0.5-2.5mm to 1.0-2.5mm

User request: zoom top/lateral in further. `Z_VIEW_MIN` changed from
`0.5e-3` to `1.0e-3` (`Z_VIEW_MAX` unchanged at `2.5e-3`). Confirmed with
the user first that this constant is shared by `render_top`, `render_lateral`,
*and* `render_xray` (unified onto one constant in an earlier round) --
they confirmed all three should narrow together rather than giving xray
its own separate window again.

Verified on all 4 demo timesteps: axes now correctly run 1.0-2.5mm on
top/lateral/xray. t1 (the earliest sample) now shows much less content
than before, since the laser had only just passed z=1.0mm at that
timestep -- expected given the narrower window, not a bug (transverse and
the stacked composite still assemble correctly around it).

### Extend laser rays in render_xray up to the domain's top boundary

User-reported (screenshot, red box on `demo_t1_xray.png`, right above the
keyhole): the near-vertical bundle of incoming rays visibly stopped
partway down the frame instead of entering from above, as if the beam
originated mid-air.

Investigated by directly parsing a `rays_laser0_<time>.vtk` file's raw
points: every one of the 360 discrete sub-rays' *recorded* path starts at
exactly y=0.200mm (checked across 10 sample rayIndex values spanning the
whole 0-359 range) -- a fixed launch plane the solver's ray-tracing model
uses internally, not the domain's actual top boundary (`ymin`, y=0 for this
case) and not the true, dynamic metal surface either. Since y=0.2mm sits
inside this view's visible y-window ([Y_DEPTH_MIN, XRAY_Y_DEPTH_MAX] =
[0.05, 0.45]mm), the recorded data alone was never going to reach the top
of the frame.

**Fix**: `_load_laser_rays()` now also returns the `rayIndex` point-data
array (identifies which discrete sub-ray each point belongs to). In
`render_xray`, for each ray, find its first recorded point (via
`np.unique(ray_idx, return_index=True)` -- points are stored in contiguous
per-ray blocks, so this gives the first-index-in-file for each ray) and add
one synthetic segment per ray from `(same x/z, ymin)` to that point, at the
same (launch) power -- so each ray visibly continues up to the domain's own
top boundary. matplotlib's existing axes clipping (`ax.set_xlim`/
`set_ylim`, already in place) cuts it off at this view's own crop exactly
like every other overlay, no separate visibility logic needed.

Verified on the exact flagged frame (t=9.744e-5s): cropped/upscaled a
side-by-side of the old vs. new render at the boxed region -- the ray
bundle above the keyhole now runs unbroken to the top edge of the frame
instead of stopping partway down. Regenerated the full 4-timestep demo set;
spot-checked t4 (a later timestep with a very different ray geometry) too.

### Fix the render_xray "black square" artifact at the X_WIDTH boundary

User confirmed the artifact with a screenshot (red box drawn on
`results/demo_t1_xray.png`, right at the bumpy top surface/plume region near
z=0.65-1.15mm) and clarified the intended fix scope via a planning
back-and-forth: keep `X_WIDTH`/`x0`/`x1` exactly as they are (the displayed
attenuation should still represent ~0.3mm of material, matching the real
sample) -- only make the *phase determination at the boundary* more robust.

Root cause confirmed: `start_state()` used `vtkPointLocator.FindClosestPoint`
-- nearest mesh *vertex* in the whole domain, not the field value actually
at `(x0, y0, z0)`. Right where the real interface passes near the fixed
`x0`/`x1` slit edge (a spatter droplet or raised bump -- exactly the boxed
region), the local mesh can be coarse enough that the nearest vertex
legitimately belongs to the wrong side of the interface, flipping that
ray's entire integrated attenuation -- visible as a blocky, cell-sized
misclassification. Deep inside uniform solid (most of the image, away from
any interface), this lookup was fine, which is why the earlier reverted fix
attempt (crossing-parity from an assumed-gas far boundary, see the entry
below) broke a *different*, much more common case instead of this one.

**Fix**: replaced the nearest-vertex lookup with proper cell interpolation
-- find the mesh cell that actually contains `(x0, y0, z0)` via a new
`vtkStaticCellLocator` (`start_cell_locator`, same class already used for
`gm_tree`/`ls_tree` elsewhere in this function) and blend `FIELD_GM`/
`FIELD_LS` using that cell's own shape-function weights
(`cell.EvaluatePosition(...)`), rather than snapping to whichever single
vertex is nearest. This is local to `(x0, y0, z0)` -- no assumption about
what lies far away -- so it carries none of the previous attempt's
regression risk. Kept the old `vtkPointLocator` purely as a fallback for
the rare case `FindCell` misses (happens ~100 times per frame out of
~370k sub-rays, always at a handful of recurring y-values -- likely
OpenFOAM's general polyhedral cells occasionally tripping up the
point-in-cell test, or a genuine coincidental alignment with an AMR
transition boundary repeating at regular z-intervals); a ray that hits the
fallback is exactly as before, not improved but not regressed either.
Logging changed from one line per fallback occurrence (too noisy for batch
runs across many frames) to a single per-frame summary count.

**Verified**:
- The exact flagged frame (t=9.744e-5s): cropped/upscaled a side-by-side of
  the old vs. new render at the boxed region -- the blocky dashed dark
  artifact along the surface bump is completely gone, replaced by a smooth
  gradient matching the surrounding attenuation.
- Regression check against the previous failed attempt's specific failure
  mode: re-rendered t=1.734e-4s and confirmed the contrast-stretch `vmin`
  stayed at exactly 0.2231 (matching the known-good pre-any-fix value, not
  the ~0.42 the reverted attempt produced) and gm crossing counts matched
  exactly (174619), confirming the solid region ahead of the laser still
  renders correctly dark.
- Regenerated the full 4-timestep demo set end to end; spot-checked t3 too
  -- no new artifacts anywhere in the image.

### Dead end: crossing-parity start-state fix, tried before the cell-interpolation fix above

Before landing on the cell-interpolation fix above, a first attempt at the
same artifact was implemented, tested, and reverted -- kept here for
history (same convention as the blue-line revert entry below, unrelated).

Implemented: extend each ray's `crossings()` cast from `xmin - margin`
(assumed gas, outside the mesh) through to `x1`, and derive the phase at
`x0` from the *parity* of gm/ls crossings before `x0`, instead of
`start_state()`'s nearest-point lookup. Rendering testrun64 at t=1.734e-4s
with this in place changed the contrast-stretch `vmin` from 0.2231 to
0.4248 and visually turned the *entire* solid substrate region ahead of the
laser (z > ~1.6mm at this timestep) pure white (0 attenuation) instead of
the correct dark solid-metal reading.

Root cause of *that*: directly sampled `alpha_smoothed` along y at
x=0,z=2.0mm (ahead of the laser) and confirmed solid metal genuinely exists
there (y from ~0.20mm to the domain bottom) -- then traced a single ray at
(y0=0.3mm, z0=2.0mm) with the new logic and found **zero** gm crossings
across the *entire* domain width, even though the far "known-gas" endpoint
and the actual solid region are unambiguously different phases. Conclusion:
the base substrate/bead spans the **full x domain width edge-to-edge** away
from the actively-melted keyhole (only the y/z-direction geometry has a
real gas/metal interface there) -- so there's no marching-cubes surface
crossing near the domain's x boundary to begin with, and assuming "just
outside the mesh = gas" silently breaks the whole parity chain for every
such ray. This is the common case *away* from the melt track, i.e. most of
the image -- so this fix, while solving the originally-hypothesized
localized coarse-cell issue in principle, was wrong far more often than the
bug it targeted. Reverted in full (`start_state`/`point_locator` restored
exactly, `_MAX_HITS_PER_RAY` back to 64) and verified byte-for-byte
equivalent crossing counts/contrast stretch to before the attempt, before
moving on to the cell-interpolation fix that actually worked (entry above).

### Colorbar title height alignment, transverse's z-minus-laser field, xray z-window + extended melt line

Four more changes to `render_view.py`, from the round of feedback right
after the previous entry below:

1. **Colorbar title vertically centered on the bar strip itself, not the
   whole label+bar image.** The title-to-the-left change (previous entry)
   centered the title against `bar_img`'s full height, which also includes
   the tick labels sitting above the strip in ParaView's own rendering --
   so the title ended up noticeably higher than the strip it labels. Added
   `_find_bar_row_range()`, which finds the strip's own row range within
   `bar_img` (the rows where nearly every pixel across the row is opaque --
   the strip is a full-width solid fill, unlike the sparse label text rows
   above it), and centers the title on that sub-range instead.
2. **`render_transverse`'s color field redefined as `z - z_laser`** (negative
   throughout the cropped region, since z <= z_cut < laser_z there) instead
   of the old `dist_behind_laser = laser_z - z` (positive). Same physical
   quantity, sign flipped: `Calculator.Function` changed from
   `(laser_z-coordsZ)*1e3` to `(coordsZ-laser_z)*1e3`, array/transfer-function
   renamed `dist_behind_laser` -> `z_minus_laser` throughout, color-scale
   bounds renamed `REL_Z_COLD/HOT` -> `Z_REL_MIN/MAX` with values negated
   (0.7/1.7 -> -1.7/-0.7) and the two `ctf.RGBPoints` entries reordered to
   stay ascending-X (white at the now-more-negative `Z_REL_MIN`, blue at
   `Z_REL_MAX`) -- same visual mapping (far from laser=white, near=blue) as
   before, just under the new sign convention. Title changed to `'z -
   z_laser'` (no units suffix, matching the field name exactly, per
   explicit request).
3. **`render_xray`'s z-window changed from scan-derived to the shared fixed
   `Z_VIEW_MIN`/`Z_VIEW_MAX` (0.5-2.5mm)**, matching top/lateral. Previously
   `zmin_crop`/`zmax_crop` were computed from the laser table's own
   scan-start/scan-end plus `CROP_PADDING` -- replaced with the same
   constants top/lateral already use. `CROP_PADDING` is now unused and was
   removed; `MELT_FRONT_OFFSET` is still used (now only for the melt-
   boundary blue line's forward edge, not the crop window -- comment
   updated). At late timesteps the laser can now sit outside the visible
   window entirely (same pre-existing behavior top/lateral already have --
   not new).
4. ~~The melt-boundary blue line in `render_xray` extended backward to the
   full crop window instead of a fixed 0.6mm-behind-the-laser cutoff~~ --
   **reverted per explicit user request** (next entry below): restored
   `MELT_WINDOW_BEHIND_FRONT = 0.6e-3` and its masking exactly as before this
   round. Items 1-3 above stand.

Verified all views individually (including early/late-timestep edge cases
for the xray z-window/laser-outside-frame case), then regenerated the full
4-timestep demo set end to end.

### Revert: xray melt-boundary blue line back to the 0.6mm-behind-front cutoff

Item 4 just above (extending the blue line to the full crop window) was
reverted per explicit user request -- back to the original
`MELT_WINDOW_BEHIND_FRONT = 0.6e-3` behavior. Items 1-3 in that entry
(colorbar title alignment, transverse's `z_minus_laser` field, xray's fixed
z-window) were not reverted, only the blue-line extent.

### Overlay the laser ray-tracing VTKs (rays + power) on render_xray, fading by power

`render_xray` now draws the solver's own multi-reflection ray-tracing output
(`VTKs/rays_laser0_<time>.vtk` -- one polyline per discrete sub-ray, broken
into many 2-point segments each carrying a `power` point-data value that
drops as the ray is absorbed/reflects) as orange line segments overlaid on
the attenuation image, alpha-faded by each segment's own power (normalized
against that frame's own peak) so a ray visibly dims out as it loses energy
along its path, rather than being drawn at constant brightness end to end.

- New `_load_laser_rays(case_dir, time_value)`: finds the right file via
  `VTKs/rays_laser0.vtk.series` (the same JSON time->filename index ParaView
  itself uses for file series -- avoids guessing at float-formatted
  filenames), reads it with `vtk.vtkPolyDataReader()` (already-imported
  `vtk`/`vtk_to_numpy`, no new dependency), and returns points/segments/
  power arrays. Returns `None` if a case has no such VTKs at all (not every
  case enables/keeps this postProcessing output -- a normal case, not an
  error -- `render_xray` logs and skips the overlay rather than failing).
- Projected onto the view's (z, y) plane (x dropped -- same "collapse the
  ray-tracing axis" treatment the attenuation image itself already uses),
  drawn via a single `matplotlib.collections.LineCollection` (350k+
  segments in a typical frame -- individual `ax.plot()` calls would be far
  too slow) with per-segment RGBA (`orange` RGB, alpha = segment's average
  endpoint power / this frame's max power).
- Explicitly pinned the axes to the image's own extent (`ax.set_xlim`/
  `set_ylim`) right after `imshow()`, before adding the ray overlay --
  needed because the rays span a much wider 3D range than the visible crop
  (multi-reflection bounces can travel across most of the domain), and
  without this, matplotlib's autoscale could grow the figure to fit them
  instead of relying on axes clipping the way the rest of the view does.
- **`RAY_MAX_OPACITY = 0.5`** caps the brightest (freshest, highest-power)
  segments at 50% alpha rather than fully opaque -- 100% was too bright/
  dominant over the attenuation image underneath (user feedback,
  2026-08-02): `alpha *= RAY_MAX_OPACITY` after the per-segment power
  normalization.

Verified on testrun64 across all 4 demo timesteps: rays render as a
converging fan into the keyhole with visibly varying opacity (denser,
brighter clusters right at the keyhole where many reflections concentrate;
fainter individual incoming rays), no performance blowup (a few seconds
added per frame), and no crash on frames/cases without ray VTKs.

### Fix mushy-outline z-window overrun, narrow colorbar overlay, transverse colorbar, title-to-the-left

Four related fixes/changes to `render_view.py`, all from the same round of
user feedback on the previous "Composite layout" round below:

1. **`render_lateral`'s x=0 mushy-boundary outline was rendering past the
   fixed 0.5-2.5mm z-window.** Root cause was twofold: `ls_slice` (the raw
   outline geometry) was built from the *unclipped* `ls_contour` with no
   z/y restriction at all (unlike `feature`, which was already properly
   clipped) -- and separately, `FRAME_MARGIN`'s zoom-out factor propagates
   through the aspect-ratio relationship between `CameraParallelScale`
   (vertical) and viewport width/height to widen the *effective* z-range the
   camera actually shows by that same ~30% factor, so even the "fixed"
   `Z_VIEW_MIN`/`Z_VIEW_MAX` constants didn't strictly bound what could
   render. Fixed narrowly by adding an explicit box `Clip` on `ls_slice`
   (bounding it to `[Y_DEPTH_MIN,Y_DEPTH_MAX] x [z_window_min,z_window_max]`)
   rather than touching `FRAME_MARGIN`, which would risk reintroducing the
   "shrinking it backfires" regression documented under task 1 above.
   Verified via cell count (264 unclipped -> 260 clipped at the same
   timestep) and visually -- the outline now correctly terminates within
   the window instead of continuing to the frame edge.
2. **Colorbars redesigned again: much narrower, and placed *on* the image
   they belong to rather than concatenated below it.** The previous round's
   "embed colorbars" design (below) concatenated the colorbar directly below
   the view content inside the same PNG, full width. Per explicit user
   correction, `_save_colorbar()` was replaced with `_overlay_colorbar()`:
   renders the colorbar at only `width_frac=0.35` of the view's own
   (already-trimmed) width, with `TransparentBackground=1`, then
   alpha-composites it directly onto the bottom-center of the view image
   -- overlapping real content there, not appended below (doesn't grow the
   canvas).
3. **Added a colorbar to `render_transverse`**, which never had one --
   `_overlay_colorbar(ctf, 'dist behind laser (mm)', ..., custom_labels=[...])`
   with explicit 3-tick labels (rather than `custom_labels=None`, which gave
   a cluttered 8+-tick auto display), matching top/lateral's convention.
4. **Colorbar title moved to the left of the bar, not stacked above it.**
   Checked ParaView 5.8's actual `ScalarBar` properties directly (via a
   throwaway probe script) before assuming: `HorizontalTitle` turned out to
   mean "keep the title horizontal rather than rotated to match a vertical
   bar", not "place it to the left" -- toggling it produced no visible
   change. `TitleLocation` doesn't exist in this ParaView version at all.
   So there's no built-in option for this layout in 5.8. Solved by rendering
   the title as its own tightly-cropped transparent PNG via matplotlib
   (`_render_title_image()`) with the scalar bar's own `Title` left empty,
   then compositing title-then-bar left-to-right into one combined image
   (`_alpha_paste()`, also used to simplify the final view-image overlay
   itself) before overlaying that combined image onto the view -- both the
   standalone colorbar file and the overlay now show the same
   title-to-the-left layout. Incidentally also fixed a pre-existing minor
   clipping bug: the old design rendered title+labels+bar all inside one
   `cb_height=100` ParaView view, which wasn't quite tall enough and clipped
   the top of the title text; the title is no longer rendered inside that
   cramped space at all.

Verified all three affected views individually, then regenerated the full
4-timestep demo set end to end.

### Composite layout: reorder, embed colorbars, trim vertical whitespace, tighten z-window, drop xray's axis title

Five related changes to `render_view.py`/`_render_stacked_video.sh`, done
together since they all touch the same per-view render/stack pipeline:

1. **Stack order**: top, transverse, lateral, xray (was top, lateral, xray,
   transverse) -- transverse now sits between the two views it's most
   directly related to.
2. **Colorbars embedded, not a separate legend row** (superseded by the
   narrower overlay-on-image design in the entry just above -- kept here for
   history). `_save_colorbar()`
   now renders at `output_png`'s own (already-trimmed) width, with
   `TransparentBackground=1`, and concatenates the result directly below
   the view content *inside `output_png` itself* -- flush, no gap. Still
   also writes a standalone colorbar file, but nothing downstream needs it
   anymore: `_render_stacked_video.sh`'s stacking step is back to a plain
   4-input vstack (no `hstack`-a-legend-row step). Verified the alpha
   channel directly (not just visually) -- checked a pixel in the middle of
   the colorbar's own blue-to-red gradient (which passes *through* a
   near-white band at its zero point) and confirmed it's fully opaque
   (`alpha=1.0`), while true background is `alpha=0.0` -- confirms
   `TransparentBackground=1` derives alpha from "was anything drawn here",
   not a white-pixel chroma-key, which would have wrongly punched a hole
   through that gradient.
3. **Vertical whitespace trim** (new `_trim_vertical_whitespace()` for
   top/lateral, `_trim_vertical_whitespace_uniform()` for transverse's 3
   panels together, since they must stay the same height to `hstack`).
   Simpler than the horizontal trim -- no min-run-length trick needed, see
   task 1's notes above for why.
4. **Fixed top/lateral z-window: 0.5mm-2.5mm** (was derived from each
   case's own scan range via `CROP_PADDING`/`MELT_FRONT_OFFSET`, now fixed
   `Z_VIEW_MIN`/`Z_VIEW_MAX` constants). Zooms in for more readable detail,
   at the cost of generality -- hardcoded to roughly testrun64's own scan
   range, won't auto-fit a case with a substantially different one.
   `lateral_xray.py`'s own z-window is untouched (still scan-derived) --
   this was scoped to top/lateral only.
5. **`render_xray()`: dropped the `z - coord (mm)` x-axis title** (tick
   labels alone still convey the scale; `bbox_inches='tight'` reclaims the
   space automatically).

Verified all 4 `--view` modes individually after the changes, then
regenerated a full 4-timestep demo set (`results/demo_t{1..4}_*.png`) end
to end through the updated stacking script.

### Merge the 4 view scripts into one, selected by a CLI flag

`lateral_screenshot.py`, `lateral_xray.py`, `top_screenshot.py`, and
`transverse_screenshot.py` merged into `results/render_view.py`, one
`render_<view>()` function each, selected via
`pvpython render_view.py --view={top,lateral,xray,transverse} <case.foam>
<time> <output.png> [<output.pvsm>]`. `--view=xray` included despite being a
fundamentally different technique (numpy ray-tracing + Beer-Lambert
attenuation, no ParaView `Show()`/`Render()` at all) -- per the original
request, one consistent entry point across all four views; it just doesn't
call any of the ParaView-render-specific shared helpers the other three do.

Shared code factored into module-level helpers/constants: `log()`,
`_load_laser_time_vs_position()`, `_laser_z_at()`,
`OFFSETS_BEHIND_LASER`/`GM_RIM_COLORS`/`GM_RIM_LINE_WIDTH`, the crop-window
constants, `_save_colorbar()`, `_trim_side_whitespace()`,
`_draw_offset_markers()` (the transverse-cut marker-line loop, parameterized
by a `make_endpoints(z_cut)` callback since top/lateral fix different axes
at the marker position). `lateral_xray.py`'s own `XRAY_Y_DEPTH_MAX = 0.45e-3`
was deliberately *not* unified with the shared `Y_DEPTH_MIN/MAX` the other
three use -- it's a genuinely different value, not a duplicate.

**Verification**: rendered all 4 views on testrun64 at t=1.974e-4s through
the merged script and md5sum-compared against the pre-merge scripts' output
-- byte-identical in every case (top, lateral, transverse, xray all matched
exactly, including their colorbar files). Also ran the updated
`_render_stacked_video.sh` (which now calls `render_view.py --view=...`
instead of the 4 separate scripts) end-to-end on 2 real timesteps.

**Found and fixed a real, pre-existing bug** while doing that end-to-end
run (unrelated to the merge itself -- present since the colorbar-legend-row
change): the stacking `ffmpeg` command split the legend row's width as
`W2 = W/2` for *both* colorbar halves, which silently only worked when `W`
was even -- for an odd `W` (2599px, hit on this run), `2*W2 = W-1`, one
pixel short of the other rows' width, and `vstack` rejected the mismatch.
Fixed by computing `W1 = W - W2` for the first half instead of reusing `W2`,
so the two halves always sum to exactly `W`. Separately, the final mp4
encode step failed on that same odd width (`libx264` requires even
width/height) -- fixed by adding `-vf scale=trunc(iw/2)*2:trunc(ih/2)*2` to
crop the rare trailing odd pixel before encoding.

The original 4 scripts were removed (`git rm`) now that `render_view.py`
verifiably replaces them -- recoverable via git history if needed
(`git log -- results/<old-script>.py`).

### Show the transverse cut lines on the lateral view too

`lateral_screenshot.py` now draws the same 3 colored vertical lines as
`top_screenshot.py` (`z = laser_z - offset`, `OFFSETS_BEHIND_LASER =
[0.7mm, 0.6mm, 0.5mm]`, red/orange/yellow via `GM_RIM_COLORS`), spanning the
full `Y_DEPTH_MIN..Y_DEPTH_MAX` height -- z is the screen-horizontal axis in
this view too (forward=-x, up=-y), so constant-z cuts render as vertical
lines here as well, same as on the top view. Coexists cleanly with the x=0
mushy-boundary outline (Done just below) already in this script.

Reused the same `x = xmax + margin` "nearest the camera" placement already
established for that outline's visibility fix (refactored into a shared
`x_marker` variable) -- confirmed the occlusion direction empirically
matches the analysis: `lateral_screenshot.py`'s camera sits at large
positive x looking toward smaller x, so larger x is nearer, same as
predicted.

Verified on testrun64 at two timesteps (saved to
`results/lateral_markers_test_{mid,early}.png`): the computed `z_cut` values
match `top_screenshot.py`'s own markers exactly at the same case/time
(0.985/1.085/1.185mm at t=1.974e-4s; 0.385/0.485/0.585mm at the early
timestep), and the early-timestep edge case (offsets pushing some cuts
before the scan start) doesn't error, consistent with `top_screenshot.py`'s
own earlier verification.

### Add an x=0 mushy-boundary outline to lateral_screenshot.py

`lateral_screenshot.py` now draws a single black outline marking the melt-
pool (mushy-zone) cross-section at x=0 (the sample's center depth plane),
overlaid on the existing gas/metal render. Mushy-surface definition matches
`transverse_screenshot.py`/`lateral_xray.py` exactly (`FIELD_LS='T'` at
`LS_TSOLIDUS=840.0`, `metal_only` via scalar `Clip` not `Threshold`); the
rim curve itself comes from an exact-plane `Slice` at `x=0` (`Normal=[1,0,0]`),
same "exact intersection, not a rendered clip" technique as
`transverse_screenshot.py`'s `slice_at_cut()`.

Visibility: rather than explicitly clipping the outline to the crop window
(as originally planned), it turned out unnecessary -- the camera's own
viewport already restricts what's visible, same as already relied on for
`top_screenshot.py`'s marker lines. The real problem was occlusion: the x=0
slice sits *inside* the solid, hidden behind the surface's near (largest-x)
face. Fixed with the same trick `top_screenshot.py` uses for its transverse-
cut markers -- since this is an orthographic camera looking straight down x,
screen position depends only on (y, z), never x, so a `Transform` filter
translates the sliced curve from x=0 to `x = xmax + margin` (nearest the
camera). This makes it the frontmost geometry without changing its shape on
screen at all.

Used `transverse_screenshot.py`'s own styling (bold black, `LineWidth=4.0`)
as the default rather than blocking on a separate design round -- verified
by rendering (see below) that it reads clearly against the blue/red
x-colored surface.

Verified on testrun64 at two timesteps (saved to
`results/lateral_outline_test_{mid,early}.png`): both run without errors,
the outline is non-empty (264 cells at t=1.974e-4s, 224 cells at the early
timestep) and renders as a clean unoccluded curve whose shape matches
`lateral_xray.py`'s ray-traced melt-pool boundary line at the same
case/time.

### Align lateral_xray.py's gas/mushy surface logic with transverse_screenshot.py

`lateral_xray.py`'s mushy-zone (liquid/solid) surface now matches
`transverse_screenshot.py`'s construction exactly:

1. `FIELD_LS`: `'epsilon1'` -&gt; `'T'`, contoured at a new `LS_TSOLIDUS = 840.0`
   constant (K) instead of the shared `ISO_THRESHOLD = 0.5` -- `epsilon1` is
   solver-derived/clamped and contours poorly; `T` is the continuous primary
   field.
2. `metal_only`: `Threshold(Scalars=FIELD_GM, ThresholdRange=[ISO_THRESHOLD,
   1.0])` -&gt; `Clip(ClipType=None, Scalars=FIELD_GM, Value=ISO_THRESHOLD,
   Invert=0)` -- interpolated clip instead of a blocky, cell-snapped
   boundary.
3. `start_state()`'s seed check: `ls_start_arr[pid] >= ISO_THRESHOLD` -&gt;
   `>= LS_TSOLIDUS`, matching the new field/threshold.
4. Updated the `_self_test_locator()` docstring, which previously explained
   a degenerate-sliver mechanism specific to the old Threshold-based
   `metal_only` (cut boundary snapped to Cartesian mesh faces) -- noted that
   mechanism shouldn't recur with the new Clip-based construction, while
   keeping the multi-candidate self-test itself (degenerate ray/triangle
   cases can still occur in principle).

The gas/metal contour (`gm_contour`) was already identical between the two
scripts and needed no change. The v3 continuous-attenuation boundary-tracing
method (`MU_SOLID_AUX`/`MU_LIQUID_AUX`, `crossings()`) is unaffected -- it
operates on whatever surface `ls_tree` is built from, agnostic to which field
produced it.

Verified on testrun64 at two timesteps (saved to
`results/xray_aligned_test_{mid,early}.png`): both run without errors, the
locator self-test finds a crossing via cell 0 on the first try for both
`gm_tree` and `ls_tree` (no degenerate-sliver retries needed, consistent
with the Clip-based boundary no longer being mesh-face-snapped), and the
resulting melt-pool boundary line is smooth and reads consistently with
`transverse_screenshot.py`'s mushy-zone rim at the same case/time.

### Separate colorbars from the images

`lateral_screenshot.py` and `top_screenshot.py` no longer bake a scalar bar
into the main render. Each now also writes a standalone colorbar PNG,
derived from `output_png`'s filename by stripping any trailing `_t<time>`
suffix (falls back to appending `_colorbar` before the extension if no such
suffix is present, e.g. for one-off manual invocations) -- so repeated
per-frame calls in a batch (as in `_render_stacked_video.sh`) overwrite the
same one file with identical content, rather than each producing a
redundant copy (the "one colorbar file per case" decision).

Implementation notes:
- The colorbar is rendered in its own dedicated `RenderView`
  (`ViewSize=[1200, 220]`), not the main view -- `GetScalarBar(ctf, cb_view)`
  with `.Visibility = 1` set directly, since the usual
  `disp.SetScalarBarVisibility(view, True)` convenience method needs a data
  representation shown in that view, and this one has none.
- Had to explicitly set `colorbar.ComponentTitle = ''` -- without an
  associated representation to infer a scalar (no-component) array from,
  ParaView defaulted to appending "Component" to the title (e.g. "y (um)
  Component"), which didn't happen in the old embedded version.
- Font-size multiplier tuned back down to match the original embedded
  look (`* 3 * 0.5`) once the view was sized appropriately (1200px wide) --
  an earlier attempt at a smaller 600x150 view with a larger multiplier
  produced overlapping labels.
- `_render_stacked_video.sh` updated: each stacked frame now has a 5th row
  (a "legend" combining the top-view and lateral-view colorbars side by
  side via `hstack`, appended below the 4 view panels via `vstack`), so the
  color-scale context isn't lost from the composite/video even though it's
  no longer baked into the individual panel images.
- Verified end-to-end on testrun64 at t=1.974e-4s: both standalone colorbar
  images render cleanly (correct title, no overlap), the main renders are
  clean (no baked-in bar), and the updated 6-input stacking filter produces
  a correct 5-row composite.

### Tune the transverse cut-offset locations

`OFFSETS_BEHIND_LASER` changed from `[0.5mm, 0.3mm, 0.1mm]` to
`[0.7mm, 0.6mm, 0.5mm]` (red/orange/yellow respectively -- user decision,
2026-08-01) in both `transverse_screenshot.py` and `top_screenshot.py`.

Verified by rendering top + transverse views at 3 sample timesteps on
testrun64 (saved to `results/sample_{early,mid,late}_{top,transverse}.png`):

- **mid/late scan**: all three cuts land inside the melt track, mushy-zone
  (black) rim visible in all three transverse panels -- reads well.
- **early scan (t=9.744e-5s, laser only 0.585mm into the 0.5mm-2.9mm scan)**:
  the two larger offsets (0.7mm, 0.6mm) push their cut z *before* the scan's
  own start -- the red panel in particular lands in untouched plate with
  essentially no melt, the orange panel only just catches the very start of
  the keyhole. This is the same "reads empty behind the melt-affected zone"
  behavior the script's own header already documents as expected (not a
  bug), just more pronounced now for early frames specifically, since these
  offsets are larger than the old ones and closer to (or past) the ~0.5mm
  melt-affected-zone limit noted there.

### Show the transverse cut lines on the top view

`transverse_screenshot.py`'s 3 panel cuts (`z = laser_z - offset`,
`OFFSETS_BEHIND_LASER = [0.5mm, 0.3mm, 0.1mm]`, colored red/orange/yellow via
`GM_RIM_COLORS`) are now also drawn on `top_screenshot.py` as vertical lines
(constant z appears vertical in that view's screen space) at the same z
positions and colors. Verified against testrun64: the computed `z_cut`
values match `transverse_screenshot.py`'s own panel offsets exactly at the
same case/time, the lines track forward with the laser across timesteps, and
an early timestep (laser near scan start) doesn't error even when a line's z
would fall outside the crop window.

`top_screenshot.py` duplicates `_laser_z_at()` and the offset/color constants
from `transverse_screenshot.py` rather than sharing them via an import,
matching this codebase's existing convention (small helpers like
`_load_laser_time_vs_position` are already duplicated identically across all
four scripts, since each runs standalone via `pvpython` in Docker) -- comments
in both files flag that the values must be kept in sync.
