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
FIELD_COLOR   = 'x_coord'         # surface coloring -- lateral (x) position,
                                   # computed below via Calculator
ISO_THRESHOLD = 0.5               # isosurface value.
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


reader = OpenFOAMReader(FileName=foam_file)
reader.CellArrays = [FIELD_GM]
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

gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

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
disp.SetScalarBarVisibility(view, True)
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
colorbar = GetScalarBar(ctf, view)
colorbar.Title = 'x (um)'
colorbar.TitleColor = [0, 0, 0]
colorbar.LabelColor = [0, 0, 0]
colorbar.TitleBold = 1
colorbar.LabelBold = 1
colorbar.Orientation = 'Horizontal'
colorbar.TitleFontSize = round(colorbar.TitleFontSize * 3 * 0.5)
colorbar.LabelFontSize = round(colorbar.LabelFontSize * 3 * 0.5)
colorbar.WindowLocation = 'LowerCenter'
colorbar.ScalarBarLength = 0.25  # kept narrow so it doesn't force the
                                  # post-render side-trim (below) wider than
                                  # the actual track content needs
colorbar.UseCustomLabels = 1
colorbar.CustomLabels = [-100.0, 0.0, 100.0]
colorbar.AddRangeLabels = 0  # otherwise the actual data min/max get added
                              # as extra labels alongside the custom ones

# Camera: lateral view, looking down +x (from outside the sample toward the
# domain center), up = -y so atmosphere is "up" in the frame -- same sense
# as the extent flip in the X-ray scripts' imshow calls. CameraParallelScale
# is half the y-window height, times a small margin factor so the content
# doesn't fill the frame edge-to-edge -- that left no room for the (now
# horizontal, larger-font) colorbar overlay to sit without covering the
# track itself. ViewSize's aspect still matches the z/y window ratio, so
# the margin is applied uniformly, not just where the colorbar sits.
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
# line while still keeping real geometry and the colorbar/text (both many
# pixels tall).
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

if output_pvsm:
    SaveState(output_pvsm)
    log(f"Saved state: {output_pvsm}")
