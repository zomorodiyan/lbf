# results/ post-processing task list

Working list of improvements to the four post-processing scripts
(`lateral_screenshot.py`, `lateral_xray.py`, `top_screenshot.py`,
`transverse_screenshot.py`) and `_render_stacked_video.sh`. Add tasks here as
they come up; move to Done as they're addressed.

## Open

### 1. Trim empty borders / tighten crop ranges to shrink image sizes

Cut off unnecessary blank space so the output images (and the stacked
composite) are smaller in both pixel dimensions and file size.

**Investigated and ruled out two of the originally-planned sub-approaches**
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
  plus a colorbar legend row -- any leftover per-view margin compounds
  there, so tightening the individual views should shrink the
  composite/mp4 too.

### 2. Merge the 4 view scripts into one, selected by a CLI flag

`lateral_screenshot.py`, `lateral_xray.py`, `top_screenshot.py`, and
`transverse_screenshot.py` already duplicate a growing amount of code
verbatim -- `_load_laser_time_vs_position()`/`_laser_z_at()`, and (per the
Done entries) `OFFSETS_BEHIND_LASER`/`GM_RIM_COLORS`/`GM_RIM_LINE_WIDTH`,
plus shared crop-window constants (`X_LATERAL_MIN/MAX`, `Y_DEPTH_MIN/MAX`,
`MELT_FRONT_OFFSET`, `CROP_PADDING`). Merge the four into a single script
that picks its view via a flag, e.g. `pvpython all_views.py --view=top
<case.foam> <time> <output.png>`, so shared logic and constants live in
exactly one place instead of N near-identical copies that can silently drift
apart (the exact risk already flagged in a couple of the Done entries below
for the offset/color constants).

Open questions to settle before implementing:
- **Scope**: `lateral_screenshot.py`/`top_screenshot.py`/`transverse_screenshot.py`
  already share the most (same "contour + ParaView-render" technique,
  same crop-window constants). `lateral_xray.py` is a fundamentally
  different technique (numpy ray-tracing + Beer-Lambert attenuation, no
  ParaView `Show()`/`Render()` at all) -- decide whether it belongs in the
  same merged script as just another `--view` branch, or stays separate
  since it barely shares implementation, only concepts.
- **Invocation compatibility**: `_render_stacked_video.sh` currently runs
  each view as its own `docker run ... pvpython <script>.py ...` call --
  the merged form just needs an extra `--view=...` argument per call, not a
  structural change to that script.
- **Sequencing**: likely better done after task 1 settles (crop
  tightening), rather than merging
  first and then having to make the same in-flight changes across a bigger
  consolidated file mid-refactor -- but flagging here in case there's a
  reason to do it first instead.

## Done

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
