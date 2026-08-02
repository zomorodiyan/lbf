# Lateral (through-thickness) view of the melt pool -- normal ParaView
# screenshot variant.
#
# Companion to lateral_xray.py: that script
# synthesizes a fake "X-ray" attenuation image by ray-tracing through the
# sample. This script instead just renders the actual 3D surface geometry
# with ParaView's own renderer -- a real screenshot, no ray tracing/custom
# math -- from the same lateral camera angle (looking down the x/through-
# thickness axis, camera up = -y so the atmosphere/plume sits at the top of
# the frame and the build plate at the bottom), matching the framing used in
# tutorials/laserbeamFoam/vdep3.pvsm.
#
# Pipeline: same core technique as the X-ray scripts -- a marching-cubes
# isosurface of the solver's recorded, pre-smoothed alpha_smoothed field on
# the native mesh (not a voxel resample, which is what produced the sawtooth
# artifact those scripts had to fix) -- but here the surface is simply shown
# with ParaView's own renderer rather than ray-traced. Colored by x (the
# lateral/through-thickness coordinate -- the axis the camera looks straight
# down) rather than a physical field: since this is an orthographic profile
# view, depth along x is otherwise entirely invisible (it's the one
# direction that never affects screen position), so coloring by it is what
# actually reveals near/far structure -- e.g. spatter sitting in front of
# vs. behind the main track.
#
# Framing: "zoom/reset to data" (the GUI button, and its scripted equivalent
# ResetCamera()) fits the camera to the *bounding box of whatever's visible*.
# Naively showing the whole gas/metal surface doesn't work: that surface's
# bounding box genuinely spans the full domain, because the untouched
# surrounding plate/powder bed is geometrically *exactly* flat (every
# triangle's normal there has zero x-component, confirmed by direct
# measurement) -- so the box is technically correct but useless for framing,
# and zooming into a saved full-domain screenshot after the fact would just
# upscale a low-res corner of it.
#
# Fix: crop with a spatial Clip (cuts at a box boundary, so everything
# inside stays fully connected -- unlike thresholding on a derived scalar
# value, which was tried first and cut a fake gap wherever the surface
# *value* passed near the flat baseline, including mid-slope on a
# continuous, curved depression wall).
#
# The z-window spans the *entire* scan, start to finish (read from
# timeVsLaserPosition), identically for every frame -- not just the track
# built up so far -- so a batch of frames shares one consistent framing
# instead of continuously zooming in as the track grows; early frames just
# show a shorter track against the same (mostly flat/unmelted) backdrop
# later frames fill in. The y-window is likewise a fixed value (see
# Y_DEPTH_MIN/MAX) rather than measured per-frame, for the same reason.
#
# (Earlier attempts intermittently produced stale/wrong output across
# timesteps in testing -- root cause turned out to be a separate bug: the
# render view's own ViewTime was never set, so Render() was showing a
# stale/default timestep independent of the per-filter
# UpdatePipeline(time=...) calls. Fixed below by setting view.ViewTime
# explicitly.)
#
# Mushy-boundary (liquid/solid) outline at x=0: a single black curve marking
# the melt-pool cross-section at the sample's center depth plane, using the
# exact same mushy-surface construction as transverse_screenshot.py
# (metal_only via scalar Clip, then Contour FIELD_LS='T' at LS_TSOLIDUS --
# not epsilon1, see that script's header for why) and the same exact-plane
# Slice technique for extracting the rim curve (sidesteps depth-buffer
# ambiguity -- see transverse_screenshot.py's slice_at_cut()). This is only
# an outline overlay, not transverse_screenshot.py's full nested surface +
# filled cap treatment.
#
# Visibility trick: the x=0 slice sits inside the solid, so it would
# normally be hidden behind the opaque gas/metal surface's near (largest-x)
# face. But this is an orthographic camera looking straight down x -- screen
# position depends only on (y, z), never on x -- so translating the sliced
# curve to x = xmax + margin (nearest the camera; see the camera comment
# below for why larger x is nearer) makes it the frontmost geometry without
# changing its shape on screen at all. Same "constant-axis nudge is
# screen-position-invariant for an orthographic camera" trick
# top_screenshot.py's transverse-cut markers use, just on x instead of y.
#
# Transverse cut markers: also draws 3 colored vertical lines (constant z,
# full Y_DEPTH_MIN..Y_DEPTH_MAX height), one per transverse_screenshot.py
# panel offset -- z = laser_z - offset, same OFFSETS_BEHIND_LASER/
# GM_RIM_COLORS values as that script and as top_screenshot.py's own version
# of this (red/orange/yellow, largest-offset-behind-laser first). z is the
# screen-horizontal axis here too (forward=-x, up=-y leaves z as the
# remaining axis), so a constant-z cut renders as a vertical line, same as
# on the top view. Uses the same x=xmax+margin front-of-camera placement as
# the mushy-boundary outline above.
#
# Run via pvpython (needs paraview.simple to read the OpenFOAM case):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/lateral_screenshot.py <case.foam> <time> <output.png> [<output.pvsm>]
#
from paraview.simple import *
from paraview import servermanager
import paraview.simple
import sys
import os
import re
import time

# ParaView auto-resets the camera to fit the visible data's *3D bounding
# sphere* on every Render() call by default -- fine for isotropic data, but
# badly wrong for this very anisotropic crop (2.5mm long, 0.35mm tall): the
# sphere's diagonal is dominated by the long z-extent, so the same
# (oversized) scale gets applied to the y-direction too, leaving the actual
# content shrunk into a small central island. Disabling this lets our own
# manually-computed CameraParallelScale (see below) actually stick.
paraview.simple._DisableFirstRenderCameraReset()

_t_start = time.time()


def log(msg):
    print(f"[{time.time() - _t_start:7.1f}s] {msg}", flush=True)


# ── Parameters ───────────────────────────────────────────────────
FIELD_GM      = 'alpha_smoothed'  # gas/metal interface source (recorded field)
FIELD_LS      = 'T'               # liquid/solid (mushy-zone) interface source
                                   # -- temperature, not epsilon1 (matches
                                   # transverse_screenshot.py; see its header)
FIELD_COLOR   = 'x_coord'         # surface coloring -- lateral (x) position,
                                   # computed below via Calculator
ISO_THRESHOLD = 0.5               # gas/metal isosurface value (alpha_smoothed).
LS_TSOLIDUS   = 840.0             # K -- mushy-zone bound, per CLAUDE.md's
                                   # physical parameters (AlSi10Mg); matches
                                   # transverse_screenshot.py's own LS_TSOLIDUS
LS_COLOR = [0.0, 0.0, 0.0]        # flat color for the mushy-zone outline --
                                   # bold black, matching
                                   # transverse_screenshot.py's equivalent rim
LS_LINE_WIDTH = 4.0                # wireframe line width (px) for the outline
                                    # -- matches transverse_screenshot.py's
                                    # LS_LINE_WIDTH
OFFSETS_BEHIND_LASER = [0.7e-3, 0.6e-3, 0.5e-3]  # meters behind the laser --
                                                   # must match
                                                   # transverse_screenshot.py's
                                                   # own OFFSETS_BEHIND_LASER
                                                   # exactly, so the marker
                                                   # lines drawn here
                                                   # correspond to that
                                                   # script's actual cut
                                                   # panels
GM_RIM_COLORS = [                    # one flat color per offset, must match
    [1.0, 0.0, 0.0],                 # transverse_screenshot.py's own
    [1.0, 0.5, 0.0],                 # GM_RIM_COLORS -- red, orange, yellow,
    [1.0, 1.0, 0.0],                 # indexed in OFFSETS_BEHIND_LASER order
]
GM_RIM_LINE_WIDTH = 3.0              # wireframe line width (px) for the
                                      # marker lines
Y_DEPTH_MIN = 0.05e-3                 # y crop -- fixed across all frames (see header)
Y_DEPTH_MAX = 0.4e-3
MELT_FRONT_OFFSET = 0.05e-3           # final laser z + this = crop window's
                                       # forward edge
CROP_PADDING = 0.03e-3                # margin added behind the scan's start
VIEW_HEIGHT_PX = 500                 # output image height; width is set to
                                      # exactly match the crop window's own
                                      # z/y aspect ratio, so there's no
                                      # letterboxing
# ─────────────────────────────────────────────────────────────────

if len(sys.argv) not in (4, 5):
    print("Usage: pvpython lateral_screenshot.py <case.foam> <time> <output.png> [<output.pvsm>]")
    sys.exit(1)

foam_file, time_value, output_png = sys.argv[1], float(sys.argv[2]), sys.argv[3]
output_pvsm = sys.argv[4] if len(sys.argv) == 5 else None


def _load_laser_time_vs_position(case_dir):
    """Parse constant/timeVsLaserPosition: a list of (t (x y z)) entries."""
    path = os.path.join(case_dir, 'constant', 'timeVsLaserPosition')
    with open(path) as f:
        content = f.read()
    entries = re.findall(
        r'\(\s*([\d.eE+-]+)\s*\(\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*\)\s*\)',
        content,
    )
    table = sorted((float(t), float(x), float(y), float(z)) for t, x, y, z in entries)
    return table


def _laser_z_at(table, t):
    """Piecewise-linear interpolation of the laser's z position at time t."""
    if t <= table[0][0]:
        return table[0][3]
    if t >= table[-1][0]:
        return table[-1][3]
    for (t0, _, _, z0), (t1, _, _, z1) in zip(table[:-1], table[1:]):
        if t0 <= t <= t1:
            frac = (t - t0) / (t1 - t0)
            return z0 + frac * (z1 - z0)
    return table[-1][3]


reader = OpenFOAMReader(FileName=foam_file)
reader.CellArrays = [FIELD_GM, FIELD_LS]
reader.Createcelltopointfiltereddata = 1
reader.UpdatePipeline(time=time_value)
log("reader loaded")

merged = MergeBlocks(Input=reader)
merged.UpdatePipeline(time=time_value)
log("blocks merged")

bounds = merged.GetDataInformation().GetBounds()
xmin, xmax, ymin, ymax, zmin, zmax = bounds
log(f"Domain bounds: x=[{xmin},{xmax}] y=[{ymin},{ymax}] z=[{zmin},{zmax}]")

laser_table = _load_laser_time_vs_position(os.path.dirname(foam_file))
scan_start_z = laser_table[0][3]
scan_end_z = laser_table[-1][3]
z_window_min = scan_start_z - CROP_PADDING
z_window_max = scan_end_z + MELT_FRONT_OFFSET + CROP_PADDING
log(f"Scan z range=[{scan_start_z*1e3:.3f},{scan_end_z*1e3:.3f}]mm; "
    f"fixed crop window: y=[{Y_DEPTH_MIN*1e3:.3f},{Y_DEPTH_MAX*1e3:.3f}]mm "
    f"z=[{z_window_min*1e3:.3f},{z_window_max*1e3:.3f}]mm")

laser_z = _laser_z_at(laser_table, time_value)
log(f"Laser z={laser_z*1e3:.3f}mm at t={time_value}")

gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

# Mushy-zone (liquid/solid) surface, restricted to inside the metal first --
# same metal_only-then-contour pattern as transverse_screenshot.py (see
# header for why Clip, not Threshold, and T, not epsilon1).
metal_only = Clip(Input=merged)
metal_only.ClipType = None
metal_only.Scalars = ['POINTS', FIELD_GM]
metal_only.Value = ISO_THRESHOLD
metal_only.Invert = 0
metal_only.UpdatePipeline(time=time_value)

ls_contour = Contour(Input=metal_only)
ls_contour.ContourBy = ['POINTS', FIELD_LS]
ls_contour.Isosurfaces = [LS_TSOLIDUS]
ls_contour.UpdatePipeline(time=time_value)
ls_poly = servermanager.Fetch(ls_contour)
log(f"Mushy-zone (liquid/solid) surface: {ls_poly.GetNumberOfCells()} cells")

# Exact plane intersection at x=0 (see header) -- not a rendered clip, so
# there's no depth-buffer ambiguity to arbitrate for this 1D curve.
ls_slice = Slice(Input=ls_contour)
ls_slice.SliceType = 'Plane'
ls_slice.SliceType.Origin = [0.0, (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0, (scan_start_z + scan_end_z) / 2.0]
ls_slice.SliceType.Normal = [1.0, 0.0, 0.0]
ls_slice.UpdatePipeline(time=time_value)
ls_slice_poly = servermanager.Fetch(ls_slice)
log(f"Mushy-zone outline at x=0: {ls_slice_poly.GetNumberOfCells()} cells "
    f"(may be empty if no melt currently straddles the sample's center depth)")

# Spatial crop (see header): a y/z bounding box, full x range, fixed across
# all frames. This only cuts geometry at the box boundary -- everything
# inside stays fully connected, unlike thresholding on a derived scalar
# value.
feature = Clip(Input=gm_contour)
feature.ClipType = 'Box'
feature.ClipType.Position = [xmin, Y_DEPTH_MIN, z_window_min]
feature.ClipType.Length = [xmax - xmin, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
feature.Invert = 1  # keep the inside of the box
feature.UpdatePipeline(time=time_value)
feature_poly = servermanager.Fetch(feature)
log(f"Cropped feature: {feature_poly.GetNumberOfCells()} cells, "
    f"bounds={feature.GetDataInformation().GetBounds()}")

xcolor = Calculator(Input=feature)
xcolor.AttributeType = 'Point Data'
xcolor.ResultArrayName = FIELD_COLOR
xcolor.Function = 'coordsX*1e6'  # meters -> micrometers, so tick labels are
                                  # plain, readable numbers (-100, 0, 100)
                                  # instead of scientific notation
xcolor.UpdatePipeline(time=time_value)

y_center = (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0
z_center = (z_window_min + z_window_max) / 2.0
x_center = (xmin + xmax) / 2.0

view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
view.ViewSize = [max(1, round(VIEW_HEIGHT_PX * (z_window_max - z_window_min) / (Y_DEPTH_MAX - Y_DEPTH_MIN))), VIEW_HEIGHT_PX]
view.ViewTime = time_value  # the view has its own time state, independent of
                             # the per-filter UpdatePipeline(time=...) calls
                             # above -- without this, Render() shows whatever
                             # the view's default/stale time is instead

disp = Show(xcolor, view)
disp.Representation = 'Surface'
ColorBy(disp, ('POINTS', FIELD_COLOR))
ctf = GetColorTransferFunction(FIELD_COLOR)
# Custom diverging map, sharply transitioning at x=0, rather than a smooth
# preset (e.g. "Cool to Warm") -- a smooth map washes both sides out to
# near-white right around zero, the one place we most want blue/red to
# already read as clearly separated. Solid blue/red saturate quickly away
# from a narrow transition band instead.
xmin_um, xmax_um = xmin * 1e6, xmax * 1e6
transition = 50.0  # width of the blue->red transition band (um) -- full
                    # saturation reached at +/-50um
ctf.RGBPoints = [
    xmin_um,     0.0, 0.0, 1.0,
    -transition, 0.0, 0.0, 1.0,
    transition,  1.0, 0.0, 0.0,
    xmax_um,     1.0, 0.0, 0.0,
]
# No colorbar overlay here -- see _save_colorbar() below, which writes it to
# its own file instead (task.md: "Separate colorbars from the images").

# Mushy-boundary outline at x=0 (see header): translate the sliced curve to
# x = xmax + margin -- nearest the camera (see camera comment below), so
# it's the frontmost geometry -- without changing its (y,z) screen position
# at all, since this orthographic camera's screen position never depends on
# x. Flat black, not colored -- same ColorArrayName=['POINTS',''] pattern
# transverse_screenshot.py uses for its own flat-colored rim outlines,
# avoiding the ColorBy(rep, None) crash noted there when there's no array to
# default to.
x_marker = xmax + 0.02 * (xmax - xmin)  # shared "front of camera" x position
                                          # for both the mushy outline below
                                          # and the transverse cut markers --
                                          # nearest the camera (see camera
                                          # comment below), guaranteed
                                          # unoccluded, invisible to screen
                                          # position either way (see header)

ls_outline = Transform(Input=ls_slice)
ls_outline.Transform = 'Transform'
ls_outline.Transform.Translate = [x_marker, 0.0, 0.0]  # the slice sits at
                                    # x=0 (Origin above), so this Translate
                                    # *is* the target x coordinate, not an
                                    # additional offset
ls_outline.UpdatePipeline(time=time_value)
ls_outline_disp = Show(ls_outline, view)
ls_outline_disp.Representation = 'Wireframe'
ls_outline_disp.ColorArrayName = ['POINTS', '']
ls_outline_disp.AmbientColor = LS_COLOR
ls_outline_disp.DiffuseColor = LS_COLOR
ls_outline_disp.LineWidth = LS_LINE_WIDTH

# Transverse cut markers (see header): one Line source per offset, at
# x = x_marker so it's nearer the camera than any real geometry and
# therefore never occluded -- same flat-color/no-scalar-bar treatment as
# ls_outline_disp above (ColorArrayName=['POINTS',''], avoiding the
# ColorBy(rep, None) crash noted in transverse_screenshot.py when there's no
# array to default to).
for _offset_idx, _offset in enumerate(OFFSETS_BEHIND_LASER):
    _z_cut = laser_z - _offset
    _marker = Line(Point1=[x_marker, Y_DEPTH_MIN, _z_cut], Point2=[x_marker, Y_DEPTH_MAX, _z_cut])
    _marker_disp = Show(_marker, view)
    _marker_disp.Representation = 'Wireframe'
    _marker_disp.ColorArrayName = ['POINTS', '']
    _marker_color = GM_RIM_COLORS[_offset_idx % len(GM_RIM_COLORS)]
    _marker_disp.AmbientColor = _marker_color
    _marker_disp.DiffuseColor = _marker_color
    _marker_disp.LineWidth = GM_RIM_LINE_WIDTH
    log(f"Transverse cut marker: offset={_offset*1e3:.2f}mm behind laser -> z={_z_cut*1e3:.3f}mm")

# Camera: lateral view, looking down +x (from outside the sample toward the
# domain center), up = -y so atmosphere is "up" in the frame -- same sense
# as the extent flip in the X-ray scripts' imshow calls. CameraParallelScale
# is half the y-window height, times a small margin factor so the content
# doesn't fill the frame edge-to-edge. ViewSize's aspect still matches the
# z/y window ratio, so the margin is applied uniformly.
#
# Tried shrinking this to 1.05 now that the colorbar is separate (task.md) --
# reverted: it backfires. Zooming in (smaller margin) makes the flat plate's
# edge-on silhouette taller in *pixels*, and the post-render side-trim below
# uses a fixed _MIN_COL_HEIGHT_PX=5 threshold to tell that silhouette apart
# from real content -- so a smaller margin pushes more of the "boring"
# background past that fixed threshold, and the side-trim ends up keeping
# *more* columns, not fewer. Measured on testrun64 at t=1.974e-4s: kept
# columns went from 2503px to 3094px (lateral) and 2443px to 3017px (top)
# when margin dropped from 1.3 to 1.05 -- both larger output files, the
# opposite of this task's goal.
FRAME_MARGIN = 1.3
view.CameraParallelProjection = 1
view.CameraViewUp = [0, -1, 0]
view.CameraFocalPoint = [x_center, y_center, z_center]
view.CameraPosition = [x_center + 2.0 * (xmax - xmin), y_center, z_center]
view.CameraParallelScale = (Y_DEPTH_MAX - Y_DEPTH_MIN) / 2.0 * FRAME_MARGIN
Render(view)

SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
log(f"Saved: {output_png}")

# Trim side whitespace: the z-window (and thus ViewSize) is fixed across an
# entire batch of frames on purpose (see header -- consistent framing so
# frames stay comparable), which means an early frame with a short track
# has a lot of blank space to the right of it. Rather than shrinking the
# camera window per-frame (which would break that shared framing), crop the
# *saved image* down to its actual content column range as a final step.
# Note this does mean cropped frames are no longer pixel-aligned/same-width
# across a batch -- this trades that away for no wasted side space.
#
# "Content" can't just mean "any non-white pixel": the untouched flat plate
# is colored (by x) same as the real track, and its edge-on silhouette is a
# real, non-white pixel row spanning the *entire* z-window regardless of how
# much track exists yet (see header on the flat-plate/zero-projected-area
# finding) -- so a naive non-white check never trims anything. Requiring a
# minimum vertical run of non-white pixels per column ignores that 1px-tall
# line while still keeping real geometry (many pixels tall).
import matplotlib.image as mpimg
import numpy as np

_img = mpimg.imread(output_png)
_non_white = np.any(_img[:, :, :3] < 0.98, axis=2)
_MIN_COL_HEIGHT_PX = 5
_content_cols = np.count_nonzero(_non_white, axis=0) >= _MIN_COL_HEIGHT_PX
if _content_cols.any():
    _cmin, _cmax = np.where(_content_cols)[0][[0, -1]]
    _pad = 15
    _cmin = max(0, _cmin - _pad)
    _cmax = min(_img.shape[1] - 1, _cmax + _pad)
    mpimg.imsave(output_png, _img[:, _cmin:_cmax + 1])
    log(f"Trimmed side whitespace: kept columns [{_cmin},{_cmax}] of {_img.shape[1]}")
else:
    log("No content found for side-trim; skipping")

# Colorbar, written to its own file instead of overlaid on the render (see
# task.md: "Separate colorbars from the images"). The color scale itself
# (ctf.RGBPoints, set above) is fixed per case, not per-frame, so this file
# is the same across every frame of a batch -- derive a per-case filename by
# stripping any "_t<time>" suffix from output_png, so repeated per-frame
# calls (as in _render_stacked_video.sh) overwrite the same file rather than
# each writing their own redundant identical copy.
_cb_match = re.match(r'^(.*?)(?:_t[\d.eE+-]+)?(\.[^.]+)$', output_png)
colorbar_png = f"{_cb_match.group(1)}_colorbar{_cb_match.group(2)}"

cb_view = CreateView('RenderView')
cb_view.OrientationAxesVisibility = 0
cb_view.Background = [1, 1, 1]
cb_view.ViewSize = [1200, 220]
colorbar = GetScalarBar(ctf, cb_view)
colorbar.Visibility = 1
colorbar.Title = 'x (um)'
colorbar.ComponentTitle = ''  # otherwise ParaView appends "Component" to
                               # the title -- harmless when the color array
                               # is bound to a shown representation (as in
                               # the old embedded version), but this
                               # standalone view has no representation to
                               # infer a scalar (no-component) array from
colorbar.TitleColor = [0, 0, 0]
colorbar.LabelColor = [0, 0, 0]
colorbar.TitleBold = 1
colorbar.LabelBold = 1
colorbar.Orientation = 'Horizontal'
colorbar.TitleFontSize = round(colorbar.TitleFontSize * 3 * 0.5)
colorbar.LabelFontSize = round(colorbar.LabelFontSize * 3 * 0.5)
colorbar.WindowLocation = 'LowerCenter'
colorbar.ScalarBarLength = 0.8
colorbar.UseCustomLabels = 1
colorbar.CustomLabels = [-100.0, 0.0, 100.0]
colorbar.AddRangeLabels = 0  # otherwise the actual data min/max get added
                              # as extra labels alongside the custom ones
Render(cb_view)
SaveScreenshot(colorbar_png, cb_view, ImageResolution=cb_view.ViewSize)
Delete(cb_view)
log(f"Saved colorbar: {colorbar_png}")

if output_pvsm:
    SaveState(output_pvsm)
    log(f"Saved state: {output_pvsm}")
