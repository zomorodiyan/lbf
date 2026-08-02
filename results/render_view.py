# Merged post-processing script: all four VDEP power-sweep views
# (top, lateral, lateral X-ray, transverse), selected via --view.
#
# Was 4 separate scripts (top_screenshot.py, lateral_screenshot.py,
# lateral_xray.py, transverse_screenshot.py) that had accumulated a growing
# amount of byte-identical duplicated code -- _load_laser_time_vs_position(),
# _laser_z_at(), the transverse-cut-marker offset/color constants
# (OFFSETS_BEHIND_LASER/GM_RIM_COLORS/GM_RIM_LINE_WIDTH), shared crop-window
# constants, the colorbar-saving block, the post-render whitespace-trim
# block -- with real risk of the copies silently drifting apart (already
# flagged in a few places before this merge). Merged into one script with
# those pieces factored into shared helpers/constants, one per-view render
# function each keeping its own view-specific logic, and a --view flag to
# select which one runs.
#
# --view=xray is a fundamentally different technique from the other three
# (numpy ray-tracing + Beer-Lambert attenuation, matplotlib output -- no
# ParaView Show()/Render() at all) and shares almost nothing at the
# implementation level with top/lateral/transverse (same "contour +
# ParaView-render" technique, most of the actual duplication). It's folded
# in here anyway for one consistent entry point across all four views, per
# the original request -- render_xray() just doesn't call any of the
# ParaView-render-specific shared helpers (_save_colorbar,
# _trim_side_whitespace, _draw_offset_markers) the other three do.
#
# Run via pvpython (needs paraview.simple to read the OpenFOAM case):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/render_view.py --view={top,lateral,xray,transverse} \
#     /workspace/<case.foam> <time> <output.png> [<output.pvsm>]
#   (<output.pvsm> is accepted but ignored for --view=xray, which never had
#   a ParaView state to save.)
#
import argparse
import os
import re
import shutil
import sys
import tempfile
import time

from paraview.simple import *
from paraview import servermanager
import paraview.simple

import matplotlib
matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

import vtk
from vtk.util.numpy_support import vtk_to_numpy

# ParaView auto-resets the camera to fit the visible data's *3D bounding
# sphere* on every Render() call by default -- fine for isotropic data, but
# badly wrong for the very anisotropic crops used here (long in z, thin in
# y): the sphere's diagonal is dominated by the long z-extent, so the same
# (oversized) scale gets applied to the y-direction too, leaving the actual
# content shrunk into a small central island. Disabling this lets each
# view's own manually-computed CameraParallelScale actually stick. Irrelevant
# to --view=xray (never calls Show()/Render()) but harmless to set regardless.
paraview.simple._DisableFirstRenderCameraReset()

_t_start = time.time()


def log(msg):
    print(f"[{time.time() - _t_start:7.1f}s] {msg}", flush=True)


# ── Shared constants (used by 2+ views; view-specific constants are defined
# inside each render_*() function instead) ─────────────────────────────────
FIELD_GM = 'alpha_smoothed'   # gas/metal interface source (recorded field --
                               # one pass of fvc::average baked in by the
                               # solver, see updateProps.H)
FIELD_LS = 'T'                # liquid/solid (mushy-zone) interface source --
                               # temperature, not epsilon1. epsilon1 is
                               # *derived* from T inside the solver (TEqn.H
                               # computes it from an enthalpy correction, then
                               # clamps it to [0,1]), so it saturates to a
                               # hard 0 or 1 almost everywhere and only
                               # genuinely varies across a very thin band --
                               # exactly the kind of near-step field that
                               # contours poorly. T is the primary,
                               # continuous field the energy equation
                               # actually solves for, with no clamping, so
                               # its solidus crossing is much
                               # better-conditioned for marching cubes.
ISO_THRESHOLD = 0.5           # gas/metal isosurface value (alpha_smoothed)
LS_TSOLIDUS = 840.0           # K -- mushy-zone bound, per CLAUDE.md's
                               # physical parameters (AlSi10Mg)
LS_COLOR = [0.0, 0.0, 0.0]    # flat black for mushy-zone rim/outline curves
LS_LINE_WIDTH = 4.0           # wireframe line width (px) for those curves
OFFSETS_BEHIND_LASER = [0.7e-3, 0.6e-3, 0.5e-3]  # meters behind the laser --
                               # transverse's cross-section cut distances,
                               # also used by top/lateral for their marker
                               # lines showing where those cuts are taken.
                               # User-chosen (2026-08-01); sits at/past the
                               # ~0.5mm melt-affected-zone edge found
                               # empirically on testrun64, so the
                               # furthest-behind cut(s) may read empty at
                               # some timesteps -- expected, not a bug.
                               # Largest offset (furthest behind) first.
GM_RIM_COLORS = [              # one flat color per offset, red/orange/yellow,
    [1.0, 0.0, 0.0],            # indexed in OFFSETS_BEHIND_LASER order
    [1.0, 0.5, 0.0],
    [1.0, 1.0, 0.0],
]
GM_RIM_LINE_WIDTH = 3.0       # wireframe line width (px) for those markers
MELT_FRONT_OFFSET = 0.05e-3   # final laser z + this = a fixed z-window's
                               # forward edge (top/lateral/xray)
CROP_PADDING = 0.03e-3        # margin added behind the scan's start
                               # (top/lateral/xray)
X_LATERAL_MIN = -0.2e-3       # x crop, narrower than the full +/-0.32mm
X_LATERAL_MAX = 0.2e-3        # domain (top/transverse)
Y_DEPTH_MIN = 0.05e-3         # y crop, fixed across all frames
Y_DEPTH_MAX = 0.4e-3           # (lateral/top/transverse -- NOT xray, which
                               # uses its own wider XRAY_Y_DEPTH_MAX, see
                               # render_xray())
VIEW_HEIGHT_PX = 500           # output image height in px (per-panel for
                               # transverse); width is set to exactly match
                               # each view's own crop-window aspect ratio,
                               # so there's no letterboxing
FRAME_MARGIN = 1.3             # small margin so content doesn't fill the
                               # frame edge-to-edge (lateral/top/transverse).
                               # NOTE: tried shrinking this to 1.05 once the
                               # colorbar became a separate file (no longer
                               # needed for headroom) -- reverted, it
                               # backfires: zooming in makes the flat
                               # plate's edge-on silhouette taller in
                               # *pixels*, and the whitespace-trim below
                               # uses a *fixed* pixel-height threshold to
                               # tell that silhouette apart from real
                               # content, so a smaller margin pushes more of
                               # the background past that threshold and the
                               # trim ends up keeping *more* columns, not
                               # fewer. Measured on testrun64 at
                               # t=1.974e-4s: kept columns went from 2503px
                               # to 3094px (lateral) and 2443px to 3017px
                               # (top) when margin dropped from 1.3 to 1.05
                               # -- both larger files. (transverse has no
                               # whitespace-trim step to begin with -- its
                               # per-panel crop window is small and fixed,
                               # not scan-spanning, so there's no
                               # accumulated blank space in the first place
                               # -- so this margin only affects its zoom
                               # level, not its file size.)


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


def _save_colorbar(ctf, title, output_png):
    """Write ctf's colorbar to its own file instead of overlaying it on the
    render (see task.md: "Separate colorbars from the images"). The color
    scale is fixed per case, not per-frame, so this file is the same across
    every frame of a batch -- derive a per-case filename by stripping any
    "_t<time>" suffix from output_png, so repeated per-frame calls (as in
    _render_stacked_video.sh) overwrite the same file rather than each
    writing their own redundant identical copy.

    Rendered in its own dedicated RenderView, not the main view --
    GetScalarBar(ctf, cb_view) with .Visibility = 1 set directly, since the
    usual disp.SetScalarBarVisibility(view, True) convenience method needs a
    data representation shown in that view, and this one has none.
    ComponentTitle must be set to '' explicitly -- without an associated
    representation to infer a scalar (no-component) array from, ParaView
    otherwise appends "Component" to the title (e.g. "y (um) Component").
    """
    m = re.match(r'^(.*?)(?:_t[\d.eE+-]+)?(\.[^.]+)$', output_png)
    colorbar_png = f"{m.group(1)}_colorbar{m.group(2)}"

    cb_view = CreateView('RenderView')
    cb_view.OrientationAxesVisibility = 0
    cb_view.Background = [1, 1, 1]
    cb_view.ViewSize = [1200, 220]
    colorbar = GetScalarBar(ctf, cb_view)
    colorbar.Visibility = 1
    colorbar.Title = title
    colorbar.ComponentTitle = ''
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
    colorbar.AddRangeLabels = 0  # otherwise the actual data min/max get
                                  # added as extra labels alongside these
    Render(cb_view)
    SaveScreenshot(colorbar_png, cb_view, ImageResolution=cb_view.ViewSize)
    Delete(cb_view)
    log(f"Saved colorbar: {colorbar_png}")


def _trim_side_whitespace(output_png):
    """Crop the *saved image* down to its actual content column range.

    The z-window is fixed across an entire batch of frames on purpose
    (consistent framing so frames stay comparable), which means an early
    frame with a short track has a lot of blank space to the right of it --
    trim that off as a final step rather than shrinking the camera window
    per-frame (which would break the shared framing). Frames end up not
    pixel-aligned/same-width across a batch -- trades that away for no
    wasted side space.

    "Content" can't just mean "any non-white pixel": the untouched flat
    plate/powder bed is colored the same as the real track, and its edge-on
    silhouette is a real, non-white pixel row spanning the *entire*
    z-window regardless of how much track exists yet -- a naive non-white
    check never trims anything. Requiring a minimum vertical run of
    non-white pixels per column ignores that 1px-tall line while still
    keeping real geometry (many pixels tall).
    """
    img = mpimg.imread(output_png)
    non_white = np.any(img[:, :, :3] < 0.98, axis=2)
    min_col_height_px = 5
    content_cols = np.count_nonzero(non_white, axis=0) >= min_col_height_px
    if content_cols.any():
        cmin, cmax = np.where(content_cols)[0][[0, -1]]
        pad = 15
        cmin = max(0, cmin - pad)
        cmax = min(img.shape[1] - 1, cmax + pad)
        mpimg.imsave(output_png, img[:, cmin:cmax + 1])
        log(f"Trimmed side whitespace: kept columns [{cmin},{cmax}] of {img.shape[1]}")
    else:
        log("No content found for side-trim; skipping")


def _draw_offset_markers(view, laser_z, make_endpoints):
    """Draw the 3 colored transverse-cut marker lines (see top/lateral
    render functions' own headers) -- one Line source per
    OFFSETS_BEHIND_LASER entry, at z = laser_z - offset, colored via
    GM_RIM_COLORS. make_endpoints(z_cut) returns the (Point1, Point2) pair
    for the Line source at that z -- callers pick which axis is "fixed at
    the front-of-camera marker position" and which spans the crop window,
    since that differs between the top and lateral cameras.

    Flat-colored, not scalar-colored: ColorArrayName=['POINTS',''], not
    ColorBy(rep, None) -- the latter crashes when there's no array to
    default to (relevant here since these Line sources carry no data).
    """
    for idx, offset in enumerate(OFFSETS_BEHIND_LASER):
        z_cut = laser_z - offset
        p1, p2 = make_endpoints(z_cut)
        marker = Line(Point1=p1, Point2=p2)
        disp = Show(marker, view)
        disp.Representation = 'Wireframe'
        disp.ColorArrayName = ['POINTS', '']
        color = GM_RIM_COLORS[idx % len(GM_RIM_COLORS)]
        disp.AmbientColor = color
        disp.DiffuseColor = color
        disp.LineWidth = GM_RIM_LINE_WIDTH
        log(f"Transverse cut marker: offset={offset*1e3:.2f}mm behind laser -> z={z_cut*1e3:.3f}mm")


# ═════════════════════════════════════════════════════════════════════════
# --view=top -- bird's-eye view of the track/spatter/powder-bed layout,
# looking straight down y (build direction). Colored by y (height relative
# to the nominal surface) since that's otherwise invisible from directly
# overhead. Also draws the 3 transverse-cut marker lines (see
# _draw_offset_markers): z is the screen-horizontal axis here (forward=+y,
# up=-x), so a constant-z cut renders as a vertical line, drawn at
# y = ymin - margin (nearer the camera than any real geometry, so never
# occluded).
# ═════════════════════════════════════════════════════════════════════════
def render_top(foam_file, time_value, output_png, output_pvsm):
    FIELD_COLOR = 'y_coord'
    SURFACE_Y = 0.2e-3  # nominal flat-plate surface height (m) -- see
                         # topoSetDict's "y surface (0.2mm)"

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

    laser_z = _laser_z_at(laser_table, time_value)
    log(f"Laser z={laser_z*1e3:.3f}mm at t={time_value}")

    gm_contour = Contour(Input=merged)
    gm_contour.ContourBy = ['POINTS', FIELD_GM]
    gm_contour.Isosurfaces = [ISO_THRESHOLD]
    gm_contour.UpdatePipeline(time=time_value)
    gm_poly = servermanager.Fetch(gm_contour)
    log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

    # Spatial crop: x/z fixed across all frames, full y -- a spatial Clip
    # only cuts geometry at the box boundary, so everything inside stays
    # fully connected (unlike thresholding on a derived scalar value).
    feature = Clip(Input=gm_contour)
    feature.ClipType = 'Box'
    feature.ClipType.Position = [X_LATERAL_MIN, ymin, z_window_min]
    feature.ClipType.Length = [X_LATERAL_MAX - X_LATERAL_MIN, ymax - ymin, z_window_max - z_window_min]
    feature.Invert = 1
    feature.UpdatePipeline(time=time_value)
    feature_poly = servermanager.Fetch(feature)
    log(f"Cropped feature: {feature_poly.GetNumberOfCells()} cells, "
        f"bounds={feature.GetDataInformation().GetBounds()}")

    ycolor = Calculator(Input=feature)
    ycolor.AttributeType = 'Point Data'
    ycolor.ResultArrayName = FIELD_COLOR
    ycolor.Function = f'(coordsY-{SURFACE_Y})*1e6'  # offset from the surface
                                                      # so 0 = "at the
                                                      # surface", positive =
                                                      # deeper/recessed
    ycolor.UpdatePipeline(time=time_value)

    x_center = (X_LATERAL_MIN + X_LATERAL_MAX) / 2.0
    y_center = (ymin + ymax) / 2.0
    z_center = (z_window_min + z_window_max) / 2.0

    view = GetActiveViewOrCreate('RenderView')
    view.OrientationAxesVisibility = 0
    view.Background = [1, 1, 1]
    view.ViewSize = [max(1, round(VIEW_HEIGHT_PX * (z_window_max - z_window_min) / (X_LATERAL_MAX - X_LATERAL_MIN))), VIEW_HEIGHT_PX]
    view.ViewTime = time_value  # the view has its own time state,
                                 # independent of the per-filter
                                 # UpdatePipeline(time=...) calls above

    disp = Show(ycolor, view)
    disp.Representation = 'Surface'
    ColorBy(disp, ('POINTS', FIELD_COLOR))
    ctf = GetColorTransferFunction(FIELD_COLOR)
    # Custom diverging map, sharply transitioning at y=surface, rather than
    # a smooth preset -- a smooth map washes both sides out to near-white
    # right around zero, the one place we most want blue/red clearly
    # separated.
    ymin_off, ymax_off = (ymin - SURFACE_Y) * 1e6, (ymax - SURFACE_Y) * 1e6
    transition = 100.0  # width of the blue->red transition band (um)
    ctf.RGBPoints = [
        ymin_off,    0.0, 0.0, 1.0,
        -transition, 0.0, 0.0, 1.0,
        transition,  1.0, 0.0, 0.0,
        ymax_off,    1.0, 0.0, 0.0,
    ]

    y_marker = ymin - 0.02 * (ymax - ymin)
    _draw_offset_markers(
        view, laser_z,
        lambda z_cut: ([X_LATERAL_MIN, y_marker, z_cut], [X_LATERAL_MAX, y_marker, z_cut]),
    )

    # Camera: top-down, looking down +y, up = -x -- with forward fixed at
    # +y (must stay top-down), a right-handed camera can't independently
    # choose both "which way is up" and "which way is right": up=+x forces
    # screen-right=-z (scan direction runs backwards, mismatching the
    # lateral view); up=-x is the one choice that gives screen-right=+z,
    # matching lateral_screenshot.py's/render_lateral's own left-to-right
    # convention.
    view.CameraParallelProjection = 1
    view.CameraViewUp = [-1, 0, 0]
    view.CameraFocalPoint = [x_center, y_center, z_center]
    view.CameraPosition = [x_center, y_center - 2.0 * (ymax - ymin), z_center]
    view.CameraParallelScale = (X_LATERAL_MAX - X_LATERAL_MIN) / 2.0 * FRAME_MARGIN
    Render(view)

    SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
    log(f"Saved: {output_png}")

    _trim_side_whitespace(output_png)
    _save_colorbar(ctf, 'y (um)', output_png)

    if output_pvsm:
        SaveState(output_pvsm)
        log(f"Saved state: {output_pvsm}")


# ═════════════════════════════════════════════════════════════════════════
# --view=lateral -- through-thickness profile, looking down x. Colored by x
# (otherwise invisible from this orthographic angle) to reveal near/far
# structure. Also draws: (1) a black outline of the mushy (liquid/solid)
# boundary at x=0 (the sample's center depth plane), and (2) the same 3
# transverse-cut marker lines as --view=top. Both are translated to
# x = xmax + margin -- nearest the camera (camera sits at large +x looking
# toward -x, so larger x is nearer) -- to sit in front of the opaque
# gas/metal surface unoccluded, without changing their (y,z) screen
# position at all (orthographic projection down x never depends on x).
# ═════════════════════════════════════════════════════════════════════════
def render_lateral(foam_file, time_value, output_png, output_pvsm):
    FIELD_COLOR = 'x_coord'

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

    # Mushy-zone surface, restricted to inside the metal first. Clip (not
    # Threshold): Threshold keeps/discards whole cells, a blocky
    # cell-resolution boundary; Clip with ClipType=None clips by the
    # Scalars/Value pair directly, interpolated like Contour, matching
    # gm_contour's own precision (Invert=0 confirmed empirically to match
    # Threshold's old [ISO_THRESHOLD, 1.0] selection).
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

    # Exact plane intersection at x=0 -- not a rendered clip, so there's no
    # depth-buffer ambiguity to arbitrate for this 1D curve (same technique
    # as render_transverse's slice_at_cut()).
    ls_slice = Slice(Input=ls_contour)
    ls_slice.SliceType = 'Plane'
    ls_slice.SliceType.Origin = [0.0, (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0, (scan_start_z + scan_end_z) / 2.0]
    ls_slice.SliceType.Normal = [1.0, 0.0, 0.0]
    ls_slice.UpdatePipeline(time=time_value)
    ls_slice_poly = servermanager.Fetch(ls_slice)
    log(f"Mushy-zone outline at x=0: {ls_slice_poly.GetNumberOfCells()} cells "
        f"(may be empty if no melt currently straddles the sample's center depth)")

    # Spatial crop: y/z fixed, full x range. See render_top's own comment
    # for why Clip over a scalar Threshold.
    feature = Clip(Input=gm_contour)
    feature.ClipType = 'Box'
    feature.ClipType.Position = [xmin, Y_DEPTH_MIN, z_window_min]
    feature.ClipType.Length = [xmax - xmin, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
    feature.Invert = 1
    feature.UpdatePipeline(time=time_value)
    feature_poly = servermanager.Fetch(feature)
    log(f"Cropped feature: {feature_poly.GetNumberOfCells()} cells, "
        f"bounds={feature.GetDataInformation().GetBounds()}")

    xcolor = Calculator(Input=feature)
    xcolor.AttributeType = 'Point Data'
    xcolor.ResultArrayName = FIELD_COLOR
    xcolor.Function = 'coordsX*1e6'  # meters -> micrometers, plain tick labels
    xcolor.UpdatePipeline(time=time_value)

    y_center = (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0
    z_center = (z_window_min + z_window_max) / 2.0
    x_center = (xmin + xmax) / 2.0

    view = GetActiveViewOrCreate('RenderView')
    view.OrientationAxesVisibility = 0
    view.Background = [1, 1, 1]
    view.ViewSize = [max(1, round(VIEW_HEIGHT_PX * (z_window_max - z_window_min) / (Y_DEPTH_MAX - Y_DEPTH_MIN))), VIEW_HEIGHT_PX]
    view.ViewTime = time_value

    disp = Show(xcolor, view)
    disp.Representation = 'Surface'
    ColorBy(disp, ('POINTS', FIELD_COLOR))
    ctf = GetColorTransferFunction(FIELD_COLOR)
    xmin_um, xmax_um = xmin * 1e6, xmax * 1e6
    transition = 50.0  # width of the blue->red transition band (um)
    ctf.RGBPoints = [
        xmin_um,     0.0, 0.0, 1.0,
        -transition, 0.0, 0.0, 1.0,
        transition,  1.0, 0.0, 0.0,
        xmax_um,     1.0, 0.0, 0.0,
    ]

    # Shared "front of camera" x position for both the mushy outline and
    # the transverse cut markers below -- see this function's header.
    x_marker = xmax + 0.02 * (xmax - xmin)

    ls_outline = Transform(Input=ls_slice)
    ls_outline.Transform = 'Transform'
    ls_outline.Transform.Translate = [x_marker, 0.0, 0.0]  # ls_slice sits at
                                        # x=0 (Origin above), so this
                                        # Translate *is* the target x
                                        # coordinate, not an additional offset
    ls_outline.UpdatePipeline(time=time_value)
    ls_outline_disp = Show(ls_outline, view)
    ls_outline_disp.Representation = 'Wireframe'
    ls_outline_disp.ColorArrayName = ['POINTS', '']
    ls_outline_disp.AmbientColor = LS_COLOR
    ls_outline_disp.DiffuseColor = LS_COLOR
    ls_outline_disp.LineWidth = LS_LINE_WIDTH

    _draw_offset_markers(
        view, laser_z,
        lambda z_cut: ([x_marker, Y_DEPTH_MIN, z_cut], [x_marker, Y_DEPTH_MAX, z_cut]),
    )

    # Camera: looking down +x, up = -y so atmosphere is "up" in the frame.
    view.CameraParallelProjection = 1
    view.CameraViewUp = [0, -1, 0]
    view.CameraFocalPoint = [x_center, y_center, z_center]
    view.CameraPosition = [x_center + 2.0 * (xmax - xmin), y_center, z_center]
    view.CameraParallelScale = (Y_DEPTH_MAX - Y_DEPTH_MIN) / 2.0 * FRAME_MARGIN
    Render(view)

    SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
    log(f"Saved: {output_png}")

    _trim_side_whitespace(output_png)
    _save_colorbar(ctf, 'x (um)', output_png)

    if output_pvsm:
        SaveState(output_pvsm)
        log(f"Saved state: {output_pvsm}")


# ═════════════════════════════════════════════════════════════════════════
# --view=transverse -- 3 side-by-side cross-section panels, each looking
# down -z at a fixed distance behind the laser (OFFSETS_BEHIND_LASER), the
# classic melt-pool cross-section shot from AM literature. Each panel: the
# gas/metal surface (colored by distance behind the laser) with the mushy
# surface nested inside it (flat gray), plus exact-plane rim outlines of
# both at the cut, plus a filled cross-section cap showing the true
# cross-sectional area. See the per-panel comments below for why each
# piece is built the way it is -- this is the most involved of the four
# views and the one lateral/xray's own mushy-surface logic was aligned to
# match.
# ═════════════════════════════════════════════════════════════════════════
def render_transverse(foam_file, time_value, output_png, output_pvsm):
    CAP_SOLID_COLOR = [0.8, 0.8, 0.8]      # flat light gray, filled solid
                                            # portion of the cross-section cap
    LS_SURFACE_COLOR = [0.75, 0.75, 0.75]  # flat lighter gray, mushy-zone
                                            # surface itself (not colored)
    PANEL_GAP_PX = 20                      # white gap between panels

    # Color-scale bounds for dist_behind_laser (mm): REL_Z_COLD is the
    # first-listed (largest/furthest-from-laser) offset, REL_Z_HOT is the
    # furthest offset + 1mm -- the whole panel set shares one consistent,
    # comparable scale.
    REL_Z_COLD = OFFSETS_BEHIND_LASER[0] * 1e3
    REL_Z_HOT = max(OFFSETS_BEHIND_LASER) * 1e3 + 1.0

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
    laser_z = _laser_z_at(laser_table, time_value)
    log(f"Laser z={laser_z*1e3:.3f}mm at t={time_value}")

    gm_contour = Contour(Input=merged)
    gm_contour.ContourBy = ['POINTS', FIELD_GM]
    gm_contour.Isosurfaces = [ISO_THRESHOLD]
    gm_contour.UpdatePipeline(time=time_value)
    gm_poly = servermanager.Fetch(gm_contour)
    log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

    # Mushy-zone surface restricted to inside the metal first, same
    # Clip-then-Contour pattern as render_lateral (see its own comment).
    # metal_only also feeds cap_full below (the filled cross-section cap),
    # whose outer/gas-facing edge should align with gm_rim -- Clip's
    # interpolated boundary keeps that alignment exact, unlike a blocky
    # Threshold.
    metal_only = Clip(Input=merged)
    metal_only.ClipType = None
    metal_only.Scalars = ['POINTS', FIELD_GM]
    metal_only.Value = ISO_THRESHOLD
    metal_only.Invert = 0
    metal_only.UpdatePipeline(time=time_value)

    ls_contour = Contour(Input=metal_only)
    ls_contour.ContourBy = ['POINTS', FIELD_LS]
    ls_contour.Isosurfaces = [LS_TSOLIDUS]  # outer edge of the melt-affected
                                              # zone -- the standard
                                              # "melt pool boundary"
    ls_contour.UpdatePipeline(time=time_value)
    ls_poly = servermanager.Fetch(ls_contour)
    log(f"Mushy-zone (liquid/solid) surface: {ls_poly.GetNumberOfCells()} cells")

    x_center = (X_LATERAL_MIN + X_LATERAL_MAX) / 2.0
    y_center = (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0
    view_w = max(1, round(VIEW_HEIGHT_PX * (X_LATERAL_MAX - X_LATERAL_MIN) / (Y_DEPTH_MAX - Y_DEPTH_MIN)))

    def box_clip(input_poly, z_window_min, z_window_max):
        """x/y fixed crop window, z in [z_window_min, z_window_max]."""
        c = Clip(Input=input_poly)
        c.ClipType = 'Box'
        c.ClipType.Position = [X_LATERAL_MIN, Y_DEPTH_MIN, z_window_min]
        c.ClipType.Length = [X_LATERAL_MAX - X_LATERAL_MIN, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
        c.Invert = 1
        c.UpdatePipeline(time=time_value)
        return c

    def slice_at_cut(input_poly, z_cut):
        """Exact plane intersection at z=z_cut, then crop to the x/y window
        -- no z window needed since the slice is already flat at z_cut.
        Not a rendered clip-and-show: this sidesteps the depth-buffer
        ambiguity problem entirely, since a Slice computes the exact 1D
        intersection curve mathematically."""
        s = Slice(Input=input_poly)
        s.SliceType = 'Plane'
        s.SliceType.Origin = [x_center, y_center, z_cut]
        s.SliceType.Normal = [0.0, 0.0, 1.0]
        s.UpdatePipeline(time=time_value)
        return box_clip(s, z_cut - 1e-9, z_cut + 1e-9)

    tmp_dir = tempfile.mkdtemp(prefix="transverse_")
    panel_pngs = []

    try:
        for panel_idx, offset in enumerate(OFFSETS_BEHIND_LASER):
            z_cut = laser_z - offset
            z_center = (zmin + z_cut) / 2.0
            log(f"Panel offset={offset*1e3:.2f}mm behind laser -> z_cut={z_cut*1e3:.3f}mm, "
                f"kept z=[{zmin*1e3:.3f},{z_cut*1e3:.3f}]mm (full part)")

            # Spatial crop: x/y fixed window, z from the full domain start
            # through the cut -- with an opaque, orthographic, straight-
            # down-z camera, a fragment only shows up if it's the
            # nearest-to-camera (largest-z) hit at its (x,y), so extending
            # the kept region further back than whatever's nearest can
            # never change what's drawn -- but it does mean the untouched
            # flat plate/powder bed renders as the solid, continuous part
            # it really is (not an empty gap or knife-edge sliver).
            feature = box_clip(gm_contour, zmin, z_cut)
            feature_poly = servermanager.Fetch(feature)
            log(f"  Cropped: {feature_poly.GetNumberOfCells()} cells, "
                f"bounds={feature.GetDataInformation().GetBounds()}")

            ls_surface = box_clip(ls_contour, zmin, z_cut)
            ls_surface_poly = servermanager.Fetch(ls_surface)
            log(f"  Cropped mushy-zone surface: {ls_surface_poly.GetNumberOfCells()} cells "
                f"(may be empty this far behind the laser)")

            ls_feature = slice_at_cut(ls_contour, z_cut)
            ls_feature_poly = servermanager.Fetch(ls_feature)
            log(f"  Mushy-zone rim at cut: {ls_feature_poly.GetNumberOfCells()} cells "
                f"(may be empty this far behind the laser)")

            gm_rim = slice_at_cut(gm_contour, z_cut)
            gm_rim_poly = servermanager.Fetch(gm_rim)
            log(f"  Gas/metal rim outline at cut: {gm_rim_poly.GetNumberOfCells()} cells")

            # dist_behind_laser (mm) for coloring the gas/metal surface --
            # laser_z minus each point's own z.
            gm_zcolor = Calculator(Input=feature)
            gm_zcolor.AttributeType = 'Point Data'
            gm_zcolor.ResultArrayName = 'dist_behind_laser'
            gm_zcolor.Function = f'({laser_z}-coordsZ)*1e3'
            gm_zcolor.UpdatePipeline(time=time_value)

            # Filled cross-section cap: slice the metal *volume* (not just
            # its outer shell) at z=z_cut, then keep only the solid
            # portion. A genuinely different object from gm_rim/ls_feature
            # above (1D curves from slicing surfaces) -- this is a 2D
            # filled area from slicing the volume. Clip by scalar (T), not
            # Threshold, for the same "interpolated not cell-snapped"
            # reason as metal_only above.
            cap_full = slice_at_cut(metal_only, z_cut)
            cap_solid = Clip(Input=cap_full)
            cap_solid.ClipType = None
            cap_solid.Scalars = ['POINTS', FIELD_LS]
            cap_solid.Value = LS_TSOLIDUS
            cap_solid.Invert = 1  # keep T < Value (solid side)
            cap_solid.UpdatePipeline(time=time_value)
            cap_solid_poly = servermanager.Fetch(cap_solid)
            log(f"  Cross-section cap (solid only): {cap_solid_poly.GetNumberOfCells()} cells")

            view = CreateView('RenderView')
            view.OrientationAxesVisibility = 0
            view.Background = [1, 1, 1]
            view.ViewSize = [view_w, VIEW_HEIGHT_PX]
            view.ViewTime = time_value

            # Shown first (flat gray, not colored) -- both fully opaque, so
            # the z-buffer decides which one wins at a given pixel
            # regardless of Show() order, but gm is drawn *after* this so
            # it wins any exact-tie/z-fighting case too (gm is always the
            # true outer boundary, mush is always nested inside it).
            ls_surface_disp = Show(ls_surface, view)
            ls_surface_disp.Representation = 'Surface'
            ls_surface_disp.ColorArrayName = ['POINTS', '']
            ls_surface_disp.AmbientColor = LS_SURFACE_COLOR
            ls_surface_disp.DiffuseColor = LS_SURFACE_COLOR

            disp = Show(gm_zcolor, view)
            disp.Representation = 'Surface'
            ColorBy(disp, ('POINTS', 'dist_behind_laser'))

            ctf = GetColorTransferFunction('dist_behind_laser')
            # Smooth blue->white gradient across the full REL_Z_COLD-HOT
            # range. White is the endpoint, not a midpoint on the way to
            # red -- ParaView's default diverging interpolation would
            # otherwise wash out through white *then* continue past it.
            ctf.RGBPoints = [
                REL_Z_COLD, 0.0, 0.0, 1.0,
                REL_Z_HOT,  1.0, 1.0, 1.0,
            ]

            # Cross-section cap's solid portion, drawn last (frontmost --
            # it sits exactly at the clip window's own front face, so this
            # only matters for exact ties, not real occlusion).
            cap_solid_disp = Show(cap_solid, view)
            cap_solid_disp.Representation = 'Surface'
            cap_solid_disp.ColorArrayName = ['POINTS', '']
            cap_solid_disp.AmbientColor = CAP_SOLID_COLOR
            cap_solid_disp.DiffuseColor = CAP_SOLID_COLOR

            # No colorbar/text label overlay -- bare render.
            ls_disp = Show(ls_feature, view)
            ls_disp.Representation = 'Wireframe'
            # Not ColorBy(ls_disp, None) -- when ls_feature is empty (no
            # melt this far behind the laser, common at early
            # timesteps/far offsets), Show() has no array to default to,
            # and ColorBy(rep, None) then crashes on that lookup. Setting
            # ColorArrayName directly is the same end state without it.
            ls_disp.ColorArrayName = ['POINTS', '']
            ls_disp.AmbientColor = LS_COLOR
            ls_disp.DiffuseColor = LS_COLOR
            ls_disp.LineWidth = LS_LINE_WIDTH

            gm_rim_disp = Show(gm_rim, view)
            gm_rim_disp.Representation = 'Wireframe'
            gm_rim_disp.ColorArrayName = ['POINTS', '']
            gm_rim_color = GM_RIM_COLORS[panel_idx % len(GM_RIM_COLORS)]
            gm_rim_disp.AmbientColor = gm_rim_color
            gm_rim_disp.DiffuseColor = gm_rim_color
            gm_rim_disp.LineWidth = GM_RIM_LINE_WIDTH

            # Camera: looking down -z (from further along the scan
            # direction back toward the cut face), up = -y.
            view.CameraParallelProjection = 1
            view.CameraViewUp = [0, -1, 0]
            view.CameraFocalPoint = [x_center, y_center, z_center]
            view.CameraPosition = [x_center, y_center, z_center + 2.0 * (zmax - zmin)]
            view.CameraParallelScale = (Y_DEPTH_MAX - Y_DEPTH_MIN) / 2.0 * FRAME_MARGIN
            Render(view)

            panel_png = os.path.join(tmp_dir, f"panel_{offset*1e3:.1f}mm.png")
            SaveScreenshot(panel_png, view, ImageResolution=view.ViewSize)
            panel_pngs.append(panel_png)
            log(f"  Saved panel: {panel_png}")

            Delete(view)

        if output_pvsm:
            SaveState(output_pvsm)
            log(f"Saved state: {output_pvsm}")

        # Stack the panels side-by-side. All panels share one fixed crop
        # window, so every panel is already the same size and fully
        # populated -- no per-panel whitespace trim needed.
        imgs = [mpimg.imread(p) for p in panel_pngs]
        gap = np.ones((imgs[0].shape[0], PANEL_GAP_PX, imgs[0].shape[2]), dtype=imgs[0].dtype)
        pieces = []
        for i, im in enumerate(imgs):
            if i > 0:
                pieces.append(gap)
            pieces.append(im)
        composite = np.concatenate(pieces, axis=1)
        mpimg.imsave(output_png, composite)
        log(f"Saved composite: {output_png}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ═════════════════════════════════════════════════════════════════════════
# --view=xray -- synthetic X-ray attenuation projection, looking down x.
# Fundamentally different technique from the other three: no ParaView
# Show()/Render() at all -- a marching-cubes isosurface of alpha_smoothed
# (and, restricted to the metal, of T at TSolidus) is extracted with
# ParaView, fetched into VTK polydata, then ray-traced by hand (Beer-Lambert
# attenuation through a per-ray phase sequence found via
# vtkStaticCellLocator) and rendered with matplotlib -- not a screenshot of
# the actual 3D geometry. See the ray-crossing/self-test machinery below for
# the vtkStaticCellLocator gotchas this works around.
#
# Draws a dotted melt-pool boundary line from a *continuous* attenuation-
# ceiling test (not a boolean presence flag -- see THRESHOLD_FRAC_MELT
# below) rather than baking the liquid/solid distinction into the displayed
# grayscale itself (MU_SOLID == MU_LIQUID for the main image): physically,
# solid and liquid AlSi10Mg are nearly the same density, and baking the
# distinction into the displayed attenuation let mesh-sliver artifacts show
# up as visible patchy/blocky noise in the grayscale itself.
# ═════════════════════════════════════════════════════════════════════════
def render_xray(foam_file, time_value, output_png):
    MU_GAS = 0.0                # attenuation coeff, gas phase (per mm)
    MU_SOLID = 5.0               # attenuation coeff, solid metal (per mm)
    MU_LIQUID = 5.0               # equal to MU_SOLID -- the melt pool
                                   # boundary is shown via the overlay line
                                   # instead of a grayscale difference
    MU_SOLID_AUX = 5.0            # solid/liquid attenuation used *only* to
    MU_LIQUID_AUX = 3.5           # locate the melt-pool boundary line --
                                   # deliberately distinct so "is this ray
                                   # 100% solid" is a meaningful, continuous
                                   # test. Never used for the displayed image.
    XRAY_Y_DEPTH_MAX = 0.45e-3    # NOTE: deeper than the shared Y_DEPTH_MAX
                                   # (0.4e-3) used by top/lateral/transverse
                                   # -- matches the 5um-refined mesh band,
                                   # see topoSetDict. Keep this view-specific,
                                   # don't unify with the shared constant.
    X_WIDTH = 0.3e-3               # ray path length (through-thickness beam
                                   # direction): the simulation domain is
                                   # 0.64mm wide, but the real experimental
                                   # sample this compares against was only
                                   # ~0.29mm -- attenuating over the full
                                   # (wider) simulation domain overstates the
                                   # path length and crushes the whole
                                   # image's contrast, so the ray is cropped
                                   # to the central X_WIDTH instead.
    NY, NZ = 120, 768             # output image resolution (rays) -- not
                                   # tied to mesh cell size, purely an
                                   # image-quality choice.

    reader = OpenFOAMReader(FileName=foam_file)
    reader.CellArrays = [FIELD_GM, FIELD_LS]
    reader.Createcelltopointfiltereddata = 1  # need point data for Contour
                                               # -- the reader's own native
                                               # cell->point averaging over
                                               # real mesh connectivity, not
                                               # a resample onto an
                                               # independent fixed-pitch
                                               # grid (which produced a
                                               # sawtooth artifact in an
                                               # earlier version)
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
    zmin_crop = scan_start_z - CROP_PADDING
    zmax_crop = scan_end_z + MELT_FRONT_OFFSET + CROP_PADDING
    log(f"Crop bounds: x=[{xmin},{xmax}] (full, uncropped) "
        f"y=[{Y_DEPTH_MIN},{XRAY_Y_DEPTH_MAX}] z=[{zmin_crop},{zmax_crop}]")

    gm_contour = Contour(Input=merged)
    gm_contour.ContourBy = ['POINTS', FIELD_GM]
    gm_contour.Isosurfaces = [ISO_THRESHOLD]
    gm_contour.UpdatePipeline(time=time_value)
    gm_poly = servermanager.Fetch(gm_contour)
    log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

    # Clip by scalar (not Threshold) -- see render_lateral's own comment;
    # same reasoning, matches transverse's construction exactly.
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
    log(f"Liquid/solid surface: {ls_poly.GetNumberOfCells()} cells")

    # NOTE: vtkModifiedBSPTree was tried here for speed but its
    # IntersectWithLine "give me every crossing" overload is unimplemented
    # in this VTK build (silently returns zero intersections).
    # vtkStaticCellLocator's multi-intersection overload has the same
    # failure; only vtkOBBTree implements it, but at ~44ms/ray that's over
    # an hour/image for a surface this size. vtkStaticCellLocator's
    # *single*-hit overload does work correctly and is built for fast
    # repeated queries against one static dataset, so "all crossings" is
    # built here by repeatedly asking for the next hit and nudging past
    # it. Verified below with an explicit self-test before trusting it.
    _NUDGE = 1e-8    # meters -- ~500x smaller than the finest 5um mesh
                     # cell, steps past a found hit without skipping a
                     # genuinely separate adjacent crossing
    _HIT_TOL = 1e-9
    _MAX_HITS_PER_RAY = 64  # generous safety cap against a pathological loop

    _crossing_totals = {'GM': 0, 'LS': 0, '_TEST': 0}

    def crossings(tree, p1, p2, tag, out):
        cur_p1 = p1
        for _ in range(_MAX_HITS_PER_RAY):
            t = vtk.mutable(0.0)
            x = [0.0, 0.0, 0.0]
            pcoords = [0.0, 0.0, 0.0]
            subId = vtk.mutable(0)
            hit = tree.IntersectWithLine(cur_p1, p2, _HIT_TOL, t, x, pcoords, subId)
            if not hit:
                break
            out.append((x[0], tag))
            _crossing_totals[tag] += 1
            cur_p1 = (x[0] + _NUDGE, cur_p1[1], cur_p1[2])
            if cur_p1[0] >= p2[0]:
                break

    def _self_test_locator(tree, poly, name, n_candidates=50):
        """Fail loudly instead of silently returning zero crossings.

        Tries several cells' centroids (not just cell 0): some cells are
        degenerate for an X-direction ray by construction, not by locator
        bug -- e.g. a triangle whose 3 vertices all share the same Y sits
        exactly in a Y-constant plane, parallel to our ray direction, a
        mathematically degenerate "ray lies in the triangle's plane" case
        any ray-triangle intersection correctly reports as no-hit. One
        such unlucky pick isn't a sign the locator is broken; only failing
        on *every* candidate is. (Any remaining slivers can still leave a
        handful of isolated pixels with a slightly-wrong liquid/solid
        split -- acceptable here since we're integrating a continuous-tone
        image, not drawing a single fragile boundary line that would
        visibly kink at every bad triangle.)
        """
        ncells = poly.GetNumberOfCells()
        n_try = min(n_candidates, ncells)
        candidate_ids = sorted({int(i * (ncells - 1) / max(n_try - 1, 1)) for i in range(n_try)})
        tried = 0
        for ci in candidate_ids:
            cell = poly.GetCell(ci)
            pts = cell.GetPoints()
            n = pts.GetNumberOfPoints()
            cy = sum(pts.GetPoint(i)[1] for i in range(n)) / n
            cz = sum(pts.GetPoint(i)[2] for i in range(n)) / n
            p1 = (xmin - abs(xmax - xmin) * 1e-4, cy, cz)
            p2 = (xmax + abs(xmax - xmin) * 1e-4, cy, cz)
            out = []
            crossings(tree, p1, p2, '_TEST', out)
            tried += 1
            if out:
                log(f"{name} self-test: {len(out)} crossing(s) found via cell {ci}/{ncells}, ok")
                return
        raise RuntimeError(
            f"{name}: self-test found zero crossings across {tried} candidate "
            f"cells -- this locator class is silently failing. Aborting "
            f"instead of producing a wrong image."
        )

    gm_tree = vtk.vtkStaticCellLocator()
    gm_tree.SetDataSet(gm_poly)
    gm_tree.BuildLocator()
    log("gm_tree built")
    _self_test_locator(gm_tree, gm_poly, "gm_tree")

    ls_tree = vtk.vtkStaticCellLocator()
    ls_tree.SetDataSet(ls_poly)
    ls_tree.BuildLocator()
    log("ls_tree built")
    _self_test_locator(ls_tree, ls_poly, "ls_tree")

    # Point locator on the native mesh, to seed each ray's start-of-line phase.
    merged_data = servermanager.Fetch(merged)
    gm_start_arr = vtk_to_numpy(merged_data.GetPointData().GetArray(FIELD_GM))
    ls_start_arr = vtk_to_numpy(merged_data.GetPointData().GetArray(FIELD_LS))
    point_locator = vtk.vtkPointLocator()
    point_locator.SetDataSet(merged_data)
    point_locator.BuildLocator()
    log("point locator built")

    def start_state(x0, y0, z0):
        pid = point_locator.FindClosestPoint((x0, y0, z0))
        metal = gm_start_arr[pid] >= ISO_THRESHOLD
        liquid = metal and (ls_start_arr[pid] >= LS_TSOLIDUS)  # T >= TSolidus
        return metal, liquid

    ys = np.linspace(Y_DEPTH_MIN, XRAY_Y_DEPTH_MAX, NY)
    zs = np.linspace(zmin_crop, zmax_crop, NZ)

    mu_integral_avg = np.empty((NZ, NY))       # main (equal-mu) -- display
    mu_integral_aux_avg = np.empty((NZ, NY))   # auxiliary -- boundary only

    x_center = (xmin + xmax) / 2.0
    x0, x1 = x_center - X_WIDTH / 2.0, x_center + X_WIDTH / 2.0
    log(f"Ray x-extent (cropped to X_WIDTH): x=[{x0},{x1}]")

    # Sub-pixel supersampling: a single ray at the exact pixel center is
    # hostage to exactly where a powder particle/spatter droplet happens to
    # sit relative to the fixed x0/x1 window. A small jittered grid of
    # sub-rays per output pixel, averaged, converges toward the true
    # local-average attenuation instead, the way a real (finite-area)
    # detector pixel would.
    N_SUB = 2  # NxN sub-ray grid per output pixel
    _SUB_FRACS = [(-0.25 + 0.5 * k / (N_SUB - 1)) if N_SUB > 1 else 0.0 for k in range(N_SUB)]
    dy_pitch = (XRAY_Y_DEPTH_MAX - Y_DEPTH_MIN) / (NY - 1)
    dz_pitch = (zmax_crop - zmin_crop) / (NZ - 1)

    def _trace_ray(y0, z0):
        metal, liquid = start_state(x0, y0, z0)
        xs = []
        crossings(gm_tree, (x0, y0, z0), (x1, y0, z0), 'GM', xs)
        crossings(ls_tree, (x0, y0, z0), (x1, y0, z0), 'LS', xs)
        xs.sort(key=lambda t: t[0])

        mu_integral = 0.0
        mu_integral_aux = 0.0
        m, l = metal, liquid
        x_prev = x0
        for x_c, tag in xs + [(x1, None)]:
            seg_len_mm = (x_c - x_prev) * 1e3
            if m:
                mu_integral += (MU_LIQUID if l else MU_SOLID) * seg_len_mm
                mu_integral_aux += (MU_LIQUID_AUX if l else MU_SOLID_AUX) * seg_len_mm
            else:
                mu_integral += MU_GAS * seg_len_mm
                mu_integral_aux += MU_GAS * seg_len_mm
            if tag == 'GM':
                m = not m
            elif tag == 'LS':
                l = not l
            x_prev = x_c
        return mu_integral, mu_integral_aux

    log(f"starting ray loop: {NY}x{NZ} pixels x {N_SUB * N_SUB} sub-rays each")
    _ray_t0 = time.time()

    for zi, z0 in enumerate(zs):
        for yi, y0 in enumerate(ys):
            mu_sum = 0.0
            mu_aux_sum = 0.0
            for dzf in _SUB_FRACS:
                zc = z0 + dzf * dz_pitch
                for dyf in _SUB_FRACS:
                    yc = y0 + dyf * dy_pitch
                    mu, mu_aux = _trace_ray(yc, zc)
                    mu_sum += mu
                    mu_aux_sum += mu_aux

            n_sub = N_SUB * N_SUB
            mu_integral_avg[zi, yi] = mu_sum / n_sub
            mu_integral_aux_avg[zi, yi] = mu_aux_sum / n_sub

        if zi == 0:
            per_row = time.time() - _ray_t0
            eta = per_row * (NZ - 1)
            log(f"row 0/{NZ} done in {per_row:.2f}s -> ETA {eta / 60:.1f} min for remaining rows")
        elif zi % 64 == 0:
            elapsed = time.time() - _ray_t0
            eta = elapsed / (zi + 1) * (NZ - zi - 1)
            log(f"row {zi}/{NZ} (ETA {eta / 60:.1f} min)")

    log("ray loop done")
    log(f"Total crossings found: GM={_crossing_totals['GM']} LS={_crossing_totals['LS']} "
        f"across {NY * NZ} rays")
    if _crossing_totals['GM'] == 0:
        raise RuntimeError(
            "Zero gas/metal crossings found across the entire ray loop -- the "
            "locator is silently failing or the surfaces/bounds are wrong. "
            "Refusing to save a bogus image."
        )

    transmission = np.exp(-mu_integral_avg)
    image = transmission.T  # (y_crop, z_crop)

    # Melt pool boundary line, from a continuous-attenuation-ceiling test
    # rather than a boolean presence flag. Needs the *auxiliary*
    # mu_integral_aux (distinct MU_SOLID_AUX/MU_LIQUID_AUX) -- its ceiling
    # (100% solid) is MU_SOLID_AUX*X_WIDTH; any liquid or gas lowers it
    # below that. THRESHOLD_FRAC_MELT relaxes the ceiling to 98%, so a ray
    # needs at least ~2% of its path to be liquid/gas before the line
    # marks it (windowed near the laser below, so this only needs to
    # reject fine mesh-sliver noise, not filter distant spurious detections).
    THRESHOLD_FRAC_MELT = 0.98
    _X_WIDTH_MM = X_WIDTH * 1e3
    _melt_ceiling = MU_SOLID_AUX * _X_WIDTH_MM

    # How far behind the melt pool's leading (front) edge to keep the blue
    # line. The front edge is derived from the laser's own known position,
    # not from the (very noisy at the strict 100% ceiling) melt detection
    # itself -- mesh-sliver noise otherwise registers as "not pure solid"
    # almost everywhere, even far ahead of where the laser has reached.
    # Scans in +z, so "front" = laser z + a small forward offset, "behind"
    # is smaller z.
    MELT_WINDOW_BEHIND_FRONT = 0.6e-3  # 0.6mm

    def _threshold_bottom(mu_avg, ceiling, y_values, frac):
        """Per z-column, deepest y where mu_avg < frac * ceiling. NaN where
        the column never falls below that mark (undisturbed)."""
        threshold = frac * ceiling
        if frac >= 0.999999:
            threshold -= 1e-6 * ceiling  # float-equality safety margin only
        below = mu_avg < threshold
        bottom = np.full(mu_avg.shape[0], np.nan)
        for zi in range(mu_avg.shape[0]):
            true_idx = np.nonzero(below[zi])[0]
            if true_idx.size:
                bottom[zi] = y_values[true_idx.max()]
        return bottom

    liquid_bottom = _threshold_bottom(mu_integral_aux_avg, _melt_ceiling, ys, THRESHOLD_FRAC_MELT)

    _laser_z = _laser_z_at(laser_table, time_value)
    _melt_front_z = _laser_z + MELT_FRONT_OFFSET
    _melt_window_min_z = _melt_front_z - MELT_WINDOW_BEHIND_FRONT
    liquid_bottom[(zs < _melt_window_min_z) | (zs > _melt_front_z)] = np.nan
    log(f"Laser z={_laser_z*1e3:.3f}mm; melt front (laser+{MELT_FRONT_OFFSET*1e3:.2f}mm)="
        f"{_melt_front_z*1e3:.3f}mm; blue line kept for "
        f"z=[{_melt_window_min_z*1e3:.3f}, {_melt_front_z*1e3:.3f}] mm")

    # Beer-Lambert transmission is long-tailed (a handful of gas/keyhole
    # pixels near 1.0 can dominate a fixed [0,1] scale). One simple rule
    # for the whole image: percentile-stretch over all of it -- with the
    # ray cropped to the physically-correct X_WIDTH, the bulk floor sits at
    # a moderate transmission, so this shows both above-surface (plume) and
    # below-surface (bulk/depression/liquid) structure without special-
    # casing any region.
    vmin, vmax = np.percentile(image, [1, 99])
    log(f"Contrast stretch: vmin={vmin:.4f} vmax={vmax:.4f} (raw range {image.min():.4f}-{image.max():.4f})")

    # True z:y physical aspect ratio (1mm of z occupies the same on-image
    # distance as 1mm of y), matching how the ParaView-rendered views are
    # scaled. Figure size computed from that same ratio so there's no
    # letterboxing.
    DPI = 150
    FIG_HEIGHT_IN = 500 / DPI
    z_range_mm = (zmax_crop - zmin_crop) * 1e3
    y_range_mm = (XRAY_Y_DEPTH_MAX - Y_DEPTH_MIN) * 1e3
    fig_width_in = FIG_HEIGHT_IN * (z_range_mm / y_range_mm)

    fig, ax = plt.subplots(figsize=(fig_width_in, FIG_HEIGHT_IN), dpi=DPI)
    extent = [zmin_crop * 1e3, zmax_crop * 1e3, XRAY_Y_DEPTH_MAX * 1e3, Y_DEPTH_MIN * 1e3]
    ax.imshow(
        image,
        cmap='gray',
        vmin=vmin, vmax=vmax,
        extent=extent,
        aspect='equal',
        interpolation='bilinear',
    )
    ax.plot(zs * 1e3, liquid_bottom * 1e3, linestyle=':', color='skyblue', linewidth=1.2,
            label='melt pool bottom (liquid/solid)', zorder=2)
    ax.axhline(0.4, color='gray', alpha=0.3, linewidth=0.8, zorder=1)  # faint
                                                                        # depth
                                                                        # reference

    from matplotlib.ticker import MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(0.5))  # round mm values only
    ax.tick_params(axis='x', labelsize=21)
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.set_xlabel('z - coord (mm)', fontsize=14)
    fig.tight_layout()
    fig.savefig(output_png, bbox_inches='tight')
    log(f"Saved: {output_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Render one of the four VDEP power-sweep post-processing views.")
    parser.add_argument('--view', required=True, choices=['top', 'lateral', 'xray', 'transverse'])
    parser.add_argument('case_foam')
    parser.add_argument('time', type=float)
    parser.add_argument('output_png')
    parser.add_argument('output_pvsm', nargs='?', default=None)
    args = parser.parse_args()

    if args.view == 'top':
        render_top(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'lateral':
        render_lateral(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'transverse':
        render_transverse(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'xray':
        if args.output_pvsm:
            log("Note: --view=xray never had a ParaView state to save; ignoring output_pvsm arg.")
        render_xray(args.case_foam, args.time, args.output_png)


if __name__ == '__main__':
    main()
