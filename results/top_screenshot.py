# Top-down view of the melt pool -- normal ParaView screenshot variant.
#
# Sibling to lateral_screenshot.py: same technique (marching-cubes isosurface
# of alpha_smoothed, rendered directly by ParaView rather than ray-traced),
# but looking straight down the y-axis (build direction) from above the
# surface instead of down the x-axis from the side -- i.e. a bird's-eye view
# of the track/spatter/powder-bed layout in the x-z (lateral vs. scan) plane.
# Colored by y (the coordinate the camera now looks straight down, so it's
# otherwise invisible) rather than a physical field, offset from the nominal
# flat-plate surface height so the diverging colormap reads as "raised above
# the surface" vs. "recessed below it" -- e.g. spatter/melt bulges vs. the
# vapor-depression keyhole, both otherwise flattened to the same silhouette
# from directly overhead.
#
# Cropping/framing follows the same reasoning as lateral_screenshot.py (see
# that file for the fuller history): the z-window spans the *entire* scan,
# start to finish, identically for every frame in a batch, so frames stay
# comparable rather than continuously zooming in as the track grows. x (now
# the screen-vertical axis) is fixed to a narrower band than the full
# +/-0.32mm domain (see X_LATERAL_MIN/MAX) for a tighter frame. y doesn't
# need cropping: it's the view axis (depth into the screen) -- opaque
# surface rendering already only shows the topmost point at each (x,z), so
# nothing needs cropping there for framing to work; it only affects
# near/far clipping.
#
# Run via pvpython (needs paraview.simple to read the OpenFOAM case):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/top_screenshot.py <case.foam> <time> <output.png> [<output.pvsm>]
#
from paraview.simple import *
from paraview import servermanager
import paraview.simple
import sys
import os
import re
import time

# See lateral_screenshot.py for why this is necessary: ParaView otherwise
# auto-resets the camera to fit the visible data's 3D bounding sphere on
# every Render() call, overriding our own manually-computed
# CameraParallelScale (below).
paraview.simple._DisableFirstRenderCameraReset()

_t_start = time.time()


def log(msg):
    print(f"[{time.time() - _t_start:7.1f}s] {msg}", flush=True)


# ── Parameters ───────────────────────────────────────────────────
FIELD_GM      = 'alpha_smoothed'  # gas/metal interface source (recorded field)
FIELD_COLOR   = 'y_coord'         # surface coloring -- height relative to the
                                   # nominal surface, computed below via Calculator
ISO_THRESHOLD = 0.5               # isosurface value.
SURFACE_Y = 0.2e-3                    # nominal flat-plate surface height (m) --
                                       # see topoSetDict's "y surface (0.2mm)"
MELT_FRONT_OFFSET = 0.05e-3           # final laser z + this = crop window's
                                       # forward edge
CROP_PADDING = 0.03e-3                # margin added behind the scan's start
X_LATERAL_MIN = -0.2e-3                # x crop -- fixed across all frames,
X_LATERAL_MAX = 0.2e-3                 # narrower than the full +/-0.32mm domain
VIEW_HEIGHT_PX = 500                 # output image height; width is set to
                                      # exactly match the crop window's own
                                      # z/x aspect ratio, so there's no
                                      # letterboxing
# ─────────────────────────────────────────────────────────────────

if len(sys.argv) not in (4, 5):
    print("Usage: pvpython top_screenshot.py <case.foam> <time> <output.png> [<output.pvsm>]")
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
    f"fixed crop window: x=[{X_LATERAL_MIN*1e3:.3f},{X_LATERAL_MAX*1e3:.3f}]mm "
    f"z=[{z_window_min*1e3:.3f},{z_window_max*1e3:.3f}]mm (y uncropped)")

gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

# Spatial crop: x (lateral) and z (scan) directions, fixed across all
# frames -- see header. A spatial Clip only cuts geometry at the box
# boundary, so everything inside stays fully connected (matches
# lateral_screenshot.py's reasoning for using Clip over a scalar Threshold).
feature = Clip(Input=gm_contour)
feature.ClipType = 'Box'
feature.ClipType.Position = [X_LATERAL_MIN, ymin, z_window_min]
feature.ClipType.Length = [X_LATERAL_MAX - X_LATERAL_MIN, ymax - ymin, z_window_max - z_window_min]
feature.Invert = 1  # keep the inside of the box
feature.UpdatePipeline(time=time_value)
feature_poly = servermanager.Fetch(feature)
log(f"Cropped feature: {feature_poly.GetNumberOfCells()} cells, "
    f"bounds={feature.GetDataInformation().GetBounds()}")

ycolor = Calculator(Input=feature)
ycolor.AttributeType = 'Point Data'
ycolor.ResultArrayName = FIELD_COLOR
ycolor.Function = f'(coordsY-{SURFACE_Y})*1e6'  # meters -> micrometers,
                                                  # offset from the nominal
                                                  # surface so 0 means "at
                                                  # the surface" -- positive
                                                  # is deeper/recessed
                                                  # (matches this codebase's
                                                  # y-as-depth convention),
                                                  # negative is raised above it
ycolor.UpdatePipeline(time=time_value)

x_center = (X_LATERAL_MIN + X_LATERAL_MAX) / 2.0
y_center = (ymin + ymax) / 2.0
z_center = (z_window_min + z_window_max) / 2.0

view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
view.ViewSize = [max(1, round(VIEW_HEIGHT_PX * (z_window_max - z_window_min) / (X_LATERAL_MAX - X_LATERAL_MIN))), VIEW_HEIGHT_PX]
view.ViewTime = time_value  # the view has its own time state, independent of
                             # the per-filter UpdatePipeline(time=...) calls
                             # above -- without this, Render() shows whatever
                             # the view's default/stale time is instead

disp = Show(ycolor, view)
disp.Representation = 'Surface'
ColorBy(disp, ('POINTS', FIELD_COLOR))
disp.SetScalarBarVisibility(view, True)
ctf = GetColorTransferFunction(FIELD_COLOR)
# Custom diverging map, sharply transitioning at y=surface (0 after the
# offset above), rather than a smooth preset -- see lateral_screenshot.py
# for why: a smooth map washes both sides out to near-white right around
# zero, the one place we most want blue/red already clearly separated.
ymin_off, ymax_off = (ymin - SURFACE_Y) * 1e6, (ymax - SURFACE_Y) * 1e6
transition = 100.0  # width of the blue->red transition band (um) -- full
                     # saturation reached at +/-100um
ctf.RGBPoints = [
    ymin_off,    0.0, 0.0, 1.0,
    -transition, 0.0, 0.0, 1.0,
    transition,  1.0, 0.0, 0.0,
    ymax_off,    1.0, 0.0, 0.0,
]
colorbar = GetScalarBar(ctf, view)
colorbar.Title = 'y (um)'
colorbar.TitleColor = [0, 0, 0]
colorbar.LabelColor = [0, 0, 0]
colorbar.TitleBold = 1
colorbar.LabelBold = 1
colorbar.Orientation = 'Horizontal'
colorbar.TitleFontSize = round(colorbar.TitleFontSize * 3 * 0.5 * 0.5 * 2)
colorbar.LabelFontSize = round(colorbar.LabelFontSize * 3 * 0.5 * 0.5 * 2)
colorbar.WindowLocation = 'LowerCenter'
colorbar.ScalarBarLength = 0.25  # kept narrow -- see lateral_screenshot.py
colorbar.UseCustomLabels = 1
colorbar.CustomLabels = [-100.0, 0.0, 100.0]
colorbar.AddRangeLabels = 0  # otherwise the actual data min/max get added
                              # as extra labels alongside the custom ones

# Camera: top-down view, looking down +y (from above the surface toward the
# domain center), up = -x -- with forward fixed at +y (must stay top-down),
# a right-handed camera can't independently choose both "which way is up"
# and "which way is right": up=+x forces screen-right=-z (scan direction
# runs backwards, laser appears to move right-to-left, mismatching the
# lateral view); up=-x is the one choice that gives screen-right=+z,
# matching lateral_screenshot.py's left-to-right convention. This also
# flips which way +x points on screen, but nothing else in this stack
# establishes a canonical x-orientation to conflict with.
# CameraParallelScale is half the x-extent, times a small margin factor so
# the content doesn't fill the frame edge-to-edge -- that would leave no
# room for the (horizontal) colorbar overlay to sit without covering the
# track itself. ViewSize's aspect matches the z/x window ratio, so the
# margin is applied uniformly.
FRAME_MARGIN = 1.3
view.CameraParallelProjection = 1
view.CameraViewUp = [-1, 0, 0]
view.CameraFocalPoint = [x_center, y_center, z_center]
view.CameraPosition = [x_center, y_center - 2.0 * (ymax - ymin), z_center]
view.CameraParallelScale = (X_LATERAL_MAX - X_LATERAL_MIN) / 2.0 * FRAME_MARGIN
Render(view)

SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
log(f"Saved: {output_png}")

# Trim side whitespace -- see lateral_screenshot.py for the full reasoning
# (the fixed z-window means short-track frames have blank space to the
# right; a naive non-white check doesn't work because the untouched flat
# plate is colored too and forms a real, non-white 1px-tall line spanning
# the entire window, so a minimum-vertical-extent filter is needed to
# ignore it while still keeping real geometry and the colorbar/text).
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
