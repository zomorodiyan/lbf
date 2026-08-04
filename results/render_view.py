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
# ParaView-render-specific shared helpers (_overlay_colorbar,
# _draw_cross_section_markers_top/_draw_cross_section_frame_lateral) the
# other three do.
#
# Colorbars (top/lateral/transverse -- xray never had one) are rendered
# narrow, with a transparent background, and alpha-composited directly onto
# the view's own output PNG near the bottom-right -- not appended below
# (which would grow the canvas) and not a separate legend file included via
# some other composition step. Title is rendered separately and placed to
# the left of the bar (ParaView 5.8's own scalar bar has no such layout
# option -- see _render_title_image()'s docstring). See _overlay_colorbar()'s
# own docstring for the rest.
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
import json
import os
import re
import shutil
import subprocess
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
from matplotlib.transforms import offset_copy
import numpy as np

import vtk
from vtk.util.numpy_support import vtk_to_numpy
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- import side effect
                                          # registers the '3d' projection
                                          # name with matplotlib on older
                                          # versions that don't do it
                                          # automatically
from mpl_toolkits.mplot3d import proj3d

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
LS_TSOLIDUS = 840.0           # K -- mushy-zone bounds, per CLAUDE.md's
LS_TLIQUIDUS = 867.0           # physical parameters (AlSi10Mg)
LS_COLOR = [0.0, 0.0, 0.0]    # flat black for solidus rim/outline curves
LS_LIQUIDUS_COLOR = [0.6, 0.6, 0.6]  # flat fainter gray for liquidus
                               # rim/outline curves -- same technique as
                               # LS_COLOR/solidus, just a lighter shade so
                               # the two are visually distinct at a glance
LS_LINE_WIDTH = 4.0           # wireframe line width (px) for those curves
# Cross-section geometry, shared by render_transverse (its own detailed
# rendering, fills etc.) and the top/lateral marker overlays below (just a
# line/frame showing where those 3 cuts sit) -- single source of truth so
# the two can't drift apart. Replaces the older OFFSETS_BEHIND_LASER/
# GM_RIM_COLORS distance-behind-laser marker system (cyan/green/yellow),
# which corresponded to a since-replaced transverse design and had nothing
# to do with render_transverse's current x-position-based cuts (2026-08-03
# redesign) -- top/lateral's markers had been silently stale ever since.
CROSS_SECTION_X = [-25e-6, 0.0, 25e-6]         # m, matches render_transverse's own
CROSS_SECTION_Z_WINDOW = (1.5e-3, 2.0e-3)      # m, fixed absolute z (not laser-relative)
CROSS_SECTION_Y_WINDOW = (0.225e-3, 0.35e-3)   # m
CROSS_SECTION_MARKER_COLORS = [                # dark blue/green/red, matches
    [0.051, 0.278, 0.631],                       # render_transverse's own
    [0.145, 0.392, 0.157],                       # DARK_COLORS hexes
    [0.718, 0.110, 0.110],                       # (#0d47a1/#256428/#b71c1c)
]
GM_RIM_LINE_WIDTH = 3.0       # wireframe line width (px) for those markers
MELT_FRONT_OFFSET = 0.05e-3   # laser z + this = the melt-pool-boundary blue
                               # line's forward edge in render_xray (the
                               # crop window itself is now the shared fixed
                               # Z_VIEW_MIN/MAX below, same as top/lateral)
Z_VIEW_MIN = 1.0e-3            # top/lateral/xray's z-window -- fixed, not
Z_VIEW_MAX = 2.5e-3             # derived from this case's own scan range.
                               # User-chosen (2026-08-02, narrowed from an
                               # earlier 0.5-2.5mm) to zoom in for
                               # readability. NOTE: hardcoded to roughly
                               # testrun64's own ~0.5-2.9mm scan -- won't
                               # auto-fit a differently-scanned case, revisit
                               # if this script is pointed at one.
X_LATERAL_MIN = -0.18e-3      # x crop, narrower than the full +/-0.32mm
X_LATERAL_MAX = 0.18e-3       # domain (top/transverse) -- narrowed from
                               # +/-0.2mm, user-chosen (2026-08-02), shared
                               # between top and transverse deliberately
Y_DEPTH_MIN = 0.05e-3         # y crop, fixed across all frames
Y_DEPTH_MAX = 0.4e-3           # (lateral/top/transverse -- NOT xray, which
                               # uses its own wider XRAY_Y_DEPTH_MAX, see
                               # render_xray())
VIEW_HEIGHT_PX = 500           # output image height in px (per-panel for
                               # transverse); width is set to exactly match
                               # each view's own crop-window aspect ratio,
                               # so there's no letterboxing
FRAME_MARGIN = 1.02             # small margin (render_top/render_lateral
                               # only) so content doesn't literally touch
                               # the frame edge. Was 1.3 -- a much larger
                               # margin, back when render_top/render_lateral
                               # still cropped their saved PNG down to
                               # content afterward (see those functions' own
                               # comments, 2026-08-03: that crop is gone
                               # now, so this margin is no longer trimmed
                               # away -- it shows up directly as wasted
                               # blank space in the final image, e.g.
                               # user-reported, 2026-08-03: "top and
                               # lateral are smaller than expected"). A
                               # previous attempt to shrink this (to 1.05)
                               # was reverted for an unrelated reason that
                               # no longer applies -- it interacted badly
                               # with that same now-removed crop's
                               # fixed-pixel-height content/background
                               # discriminator, not with anything about the
                               # rendering itself. 1.02 (not 1.0 exactly)
                               # keeps a hairline buffer so antialiased
                               # edge pixels don't get clipped.


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


def _load_laser_rays(case_dir, time_value):
    """Load the laser ray-tracing VTK (`VTKs/rays_laser0_<time>.vtk`, written
    by the solver's multi-reflection ray-tracing absorption model -- one
    polyline per discrete sub-ray, broken into many 2-point segments so each
    can carry its own `power` point-data value, which drops as the ray loses
    energy to absorption/reflection along its path) closest to time_value.

    Uses the VTKs directory's own `rays_laser0.vtk.series` (the same
    JSON time->filename index ParaView writes/reads for any file series) to
    find the right file rather than guessing at float-formatted filenames.

    Returns (points (N,3) float array, segments (M,2) int array of point-
    index pairs, power (N,) float array, rayIndex (N,) int array identifying
    which discrete sub-ray each point belongs to, matched .vtk path) -- or
    None if this case has no such VTKs at all (not every case enables/keeps
    this postProcessing output, so this is a normal, expected case, not an
    error).
    """
    series_path = os.path.join(case_dir, 'VTKs', 'rays_laser0.vtk.series')
    if not os.path.exists(series_path):
        return None
    with open(series_path) as f:
        series = json.load(f)
    best = min(series['files'], key=lambda e: abs(e['time'] - time_value))
    vtk_path = os.path.join(case_dir, 'VTKs', best['name'])

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_path)
    reader.ReadAllScalarsOn()
    reader.Update()
    poly = reader.GetOutput()

    points = vtk_to_numpy(poly.GetPoints().GetData())
    power_arr = poly.GetPointData().GetArray('power')
    power = vtk_to_numpy(power_arr).astype(np.float64) if power_arr is not None else None
    ray_idx_arr = poly.GetPointData().GetArray('rayIndex')
    ray_idx = vtk_to_numpy(ray_idx_arr) if ray_idx_arr is not None else None

    lines = poly.GetLines()
    id_list = vtk.vtkIdList()
    lines.InitTraversal()
    segments = []
    while lines.GetNextCell(id_list):
        n = id_list.GetNumberOfIds()
        for k in range(n - 1):  # handles the common 2-point-per-cell case
            segments.append((id_list.GetId(k), id_list.GetId(k + 1)))  # and
                                                                        # any
                                                                        # longer
                                                                        # polylines
    segments = np.array(segments, dtype=np.int64)

    return points, segments, power, ray_idx, vtk_path


def _alpha_paste(base, overlay, x0, y0):
    """Alpha-composite `overlay` (RGBA array) onto `base` (RGBA array) with
    overlay's top-left corner at (x0, y0) in base's pixel coordinates,
    clipping to base's bounds. Shared by _overlay_colorbar for pasting both
    the title text and the bar itself onto the view image."""
    bh, bw = base.shape[:2]
    oh, ow = overlay.shape[:2]
    x0c, y0c = max(0, x0), max(0, y0)
    x1, y1 = min(bw, x0 + ow), min(bh, y0 + oh)
    if x1 <= x0c or y1 <= y0c:
        return
    ov = overlay[y0c - y0:y1 - y0, x0c - x0:x1 - x0, :]
    region = base[y0c:y1, x0c:x1, :]
    a = ov[:, :, 3:4]
    region[:, :, :3] = region[:, :, :3] * (1 - a) + ov[:, :, :3] * a
    region[:, :, 3:4] = np.maximum(region[:, :, 3:4], a)
    base[y0c:y1, x0c:x1, :] = region


def _render_title_image(title, fontsize_pt=20, dpi=100):
    """Render `title` as a standalone tightly-cropped transparent-background
    RGBA array via matplotlib. Used instead of ParaView's own scalar-bar
    Title property, which in ParaView 5.8 always stacks the title *above* a
    horizontal bar (there's no built-in "title to the left" option -- checked
    the scalar bar's exposed properties directly: HorizontalTitle turned out
    to just mean "keep the title horizontal rather than rotated to match a
    vertical bar", not "place it to the left", and TitleLocation doesn't
    exist in this ParaView version). Rendering it ourselves gives full
    control over placement.
    """
    fig = plt.figure(figsize=(2, 1), dpi=dpi)
    fig.patch.set_alpha(0)
    text = fig.text(0.5, 0.5, title, fontsize=fontsize_pt, fontweight='bold',
                     color='black', ha='center', va='center')
    fig.canvas.draw()
    pad_px = 4
    bbox = text.get_window_extent()
    fig.set_size_inches((bbox.width + 2 * pad_px) / dpi,
                         (bbox.height + 2 * pad_px) / dpi)
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return buf.astype(np.float32) / 255.0


def _render_legend_image(entries, fontsize_pt=22, dpi=100, swatch_frac=0.22, line_width=5):
    """Render a small vertical legend (one colored line swatch + label per
    row) as a standalone matplotlib RGBA array, transparent background --
    same "render through a throwaway figure, read back the pixel buffer"
    technique as _render_title_image()/_resize_image(), just multi-row.
    Rendered oversized with generous margins; caller crops it to content
    via _trim_whitespace_bbox() (safe here since the transparent margin is
    still white RGB underneath, same as domain_schematic.png's own margins).
    entries: list of (color, label) tuples, one row each, top to bottom.
    Used by render_transverse to map color -> x-position for its overlaid
    cross-sections, since there's no single colored field to hang a
    continuous colorbar off of."""
    fig = plt.figure(figsize=(3.4, 0.55 * len(entries) + 0.3), dpi=dpi)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(entries))
    for i, (color, label) in enumerate(entries):
        y = len(entries) - i - 0.5
        ax.plot([0.04, 0.04 + swatch_frac], [y, y], color=color, linewidth=line_width,
                solid_capstyle='butt')
        ax.text(0.04 + swatch_frac + 0.04, y, label, fontsize=fontsize_pt, va='center',
                ha='left', color='black')
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return buf.astype(np.float32) / 255.0


def _find_bar_row_range(bar_img, opaque_frac_threshold=0.9):
    """Return (row_start, row_end) of the actual solid color-bar strip
    within a rendered ParaView horizontal scalar-bar image (tick labels
    below the bar strip -- see _overlay_colorbar's own TextPosition
    setting -- Title=''). Found as the run of rows where nearly every pixel
    across the row is opaque -- the bar itself is a full-width solid strip,
    while the label rows only have sparse text pixels, so a high
    opaque-fraction threshold picks out just the strip regardless of
    whether the labels sit above or below it. Falls back to the full image
    height if no such run is found. Used to vertically align the title
    against the bar strip itself, not the whole label+bar image (see
    _overlay_colorbar)."""
    alpha = bar_img[:, :, 3]
    frac_opaque = (alpha > 0.05).mean(axis=1)
    rows = np.where(frac_opaque > opaque_frac_threshold)[0]
    if len(rows) == 0:
        return 0, bar_img.shape[0]
    return int(rows.min()), int(rows.max()) + 1


def _overlay_colorbar(ctf, title, output_png, custom_labels=None,
                       width_frac=0.35, cb_height=None, vert_margin=15,
                       side_margin=15, vert_margin_frac=None,
                       down_shift_chars=0, side='right', vert='top'):
    """Render ctf's colorbar *narrower* than the view (width_frac of its
    own already-trimmed width -- deliberately small, not spanning the
    image) and alpha-composite it onto the view image itself, near the
    top-right by default (user request, 2026-08-03; went bottom-right ->
    bottom-left -> top-right over the course of that same session -- pass
    side='left'/vert='bottom' for older placements), overlapping the
    actual content there -- not appended
    below (which would grow the canvas) and not a separate legend file
    included via some other composition step. Call *after* SaveScreenshot
    (and, for render_transverse, after the schematic prepend) so
    output_png's width is already final.

    cb_height default (None) is computed as a fraction of output_png's own
    width (vw) rather than a fixed pixel count, and the bar/title font
    sizes scale with it -- top/lateral/transverse/xray all end up at
    different native vw (their own aspect ratios and content differ), but
    _render_stacked_video.sh later scales every view's saved PNG to one
    shared width for the stacked composite. A fixed absolute cb_height
    would end up a different *effective* size in that shared-width output
    depending on each view's own native vw, making the colorbars
    inconsistent across stacked rows (user report, 2026-08-02). Anchoring
    cb_height (and its fonts) to vw instead means the eventual post-scale
    size only depends on the shared target width, not each view's own
    resolution. REF_VW/REF_CB_HEIGHT is simply top view's own historical
    width and this constant's original fixed-100px tuning, kept as the
    reference point so top's own appearance is unchanged by this switch to
    proportional sizing.

    custom_labels: fixed tick values (e.g. [-100, 0, 100] for top/lateral's
    +/-um fields); leave as None to let ParaView auto-pick labels from the
    data range instead -- needed for render_transverse's z_minus_laser,
    whose range (Z_REL_MIN..Z_REL_MAX, a few mm) has nothing to do with
    the other views' +/-100um scale.

    Still also writes a standalone colorbar file (derived per-case
    filename, stripping any "_t<time>" suffix from output_png so repeated
    per-frame calls overwrite one shared file -- see task.md: "Separate
    colorbars from the images") in case it's wanted on its own.

    TransparentBackground=1 -- alpha comes from the renderer itself (was
    anything drawn at this pixel?), not a post-hoc color chroma-key. That
    distinction matters here for two reasons: this colorbar's own gradient
    can pass *through* a near-white band (e.g. top/lateral's zero point),
    which a naive "near-white = background" post-process would incorrectly
    punch a transparent hole through -- and unlike the old embed-below
    design, this alpha is now composited onto arbitrary underlying pixels
    (real geometry, not a flat canvas), so it has to be correct, not just
    visually close enough on a white background.

    Rendered in its own dedicated RenderView, not the main view --
    GetScalarBar(ctf, cb_view) with .Visibility = 1 set directly, since the
    usual disp.SetScalarBarVisibility(view, True) convenience method needs a
    data representation shown in that view, and this one has none.
    ComponentTitle must be set to '' explicitly -- without an associated
    representation to infer a scalar (no-component) array from, ParaView
    otherwise appends "Component" to the title (e.g. "y (um) Component").

    `title` is deliberately *not* set as the scalar bar's own Title property
    (left at ''): ParaView 5.8's horizontal scalar bar always stacks its
    title above the bar with no built-in "title to the left" option, so the
    title is rendered separately via _render_title_image() and composited to
    the left of the (title-less) bar via _alpha_paste() -- see that
    function's own docstring for why.

    vert_margin_frac (fraction of vh, overrides vert_margin when given):
    was used by render_transverse (no longer calls this function -- it
    builds its own in-scene labels now) instead of a fixed pixel
    vert_margin, to clear a since-removed bottom crop applied only to the
    stacked composite, by a fixed proportion of vh regardless of this
    frame's own resolution. Left available for any future caller with a
    similar need.

    down_shift_chars: nudge the whole overlaid title+bar block down the
    screen (regardless of vert='top'/'bottom') by this many characters
    (screen-space, same "0.6 * fontsize" character-height convention used
    elsewhere in this repo, e.g. plot_domain_schematic.py's axis-label
    nudges) -- positive always moves it further down; for vert='bottom'
    that also means less clearance from the bottom edge. User request,
    2026-08-02: "move the transverse colorbar image 2.5 characters down".
    """
    view_img = mpimg.imread(output_png)
    if view_img.shape[2] == 3:
        alpha = np.ones((*view_img.shape[:2], 1), dtype=view_img.dtype)
        view_img = np.concatenate([view_img, alpha], axis=2)
    vh, vw = view_img.shape[:2]
    if vert_margin_frac is not None:
        vert_margin = round(vh * vert_margin_frac)
    cb_width = max(1, round(vw * width_frac))
    REF_VW, REF_CB_HEIGHT = 1633, 100
    if cb_height is None:
        cb_height = max(1, round(vw * (REF_CB_HEIGHT / REF_VW)))
    font_scale = cb_height / REF_CB_HEIGHT

    m = re.match(r'^(.*?)(?:_t[\d.eE+-]+)?(\.[^.]+)$', output_png)
    colorbar_png = f"{m.group(1)}_colorbar{m.group(2)}"

    cb_view = CreateView('RenderView')
    cb_view.OrientationAxesVisibility = 0
    cb_view.Background = [1, 1, 1]  # irrelevant once SaveScreenshot below
                                     # renders with TransparentBackground=1
                                     # -- kept opaque here only so Render()
                                     # itself has a defined background;
                                     # never actually written to the file
    cb_view.ViewSize = [cb_width, cb_height]
    colorbar = GetScalarBar(ctf, cb_view)
    colorbar.Visibility = 1
    colorbar.Title = ''  # rendered separately and placed to the left instead
                          # -- see _render_title_image()'s docstring for why
    colorbar.ComponentTitle = ''
    colorbar.TitleColor = [0, 0, 0]
    colorbar.LabelColor = [0, 0, 0]
    colorbar.TitleBold = 1
    colorbar.LabelBold = 1
    colorbar.Orientation = 'Horizontal'
    colorbar.TitleFontSize = round(colorbar.TitleFontSize * 3 * 0.5 * font_scale)
    colorbar.LabelFontSize = round(colorbar.LabelFontSize * 3 * 0.5 * font_scale)
    label_font_size = colorbar.LabelFontSize  # captured for down_shift_chars
                                               # below -- cb_view/colorbar
                                               # get Delete()d before then
    # Tick value labels below the bar strip, not above it (the default) --
    # user request, 2026-08-02.
    colorbar.TextPosition = 'Ticks left/bottom, annotations right/top'
    colorbar.WindowLocation = 'LowerCenter'
    colorbar.ScalarBarLength = 0.8
    if custom_labels is not None:
        colorbar.UseCustomLabels = 1
        colorbar.CustomLabels = custom_labels
        colorbar.AddRangeLabels = 0  # otherwise the actual data min/max
                                      # get added as extra labels alongside
                                      # these
    Render(cb_view)
    SaveScreenshot(colorbar_png, cb_view, ImageResolution=cb_view.ViewSize,
                    TransparentBackground=1)
    Delete(cb_view)

    # Compose title (rendered standalone) to the left of the title-less bar
    # PNG just saved, then overwrite colorbar_png with that combined image --
    # so the standalone file matches what actually gets overlaid below.
    # Vertically aligned to the *bar strip itself* (via _find_bar_row_range),
    # not the whole bar_img -- bar_img also includes the tick labels sitting
    # above the strip, and centering on the whole thing would put the title
    # noticeably higher than the strip it labels.
    bar_img = mpimg.imread(colorbar_png)
    title_img = _render_title_image(title, fontsize_pt=round(colorbar.LabelFontSize * 0.9))
    gap_px = max(4, round(cb_width * 0.02))
    th, tw = title_img.shape[:2]
    bh, bw = bar_img.shape[:2]
    r0, r1 = _find_bar_row_range(bar_img)
    bar_strip_center = (r0 + r1) / 2.0

    y_bar, y_title = 0.0, bar_strip_center - th / 2.0
    top = min(y_bar, y_title)
    y_bar, y_title = y_bar - top, y_title - top
    block_h = int(np.ceil(max(y_bar + bh, y_title + th)))
    block_w = tw + gap_px + bw
    combined = np.zeros((block_h, block_w, 4), dtype=np.float32)
    _alpha_paste(combined, title_img, 0, round(y_title))
    _alpha_paste(combined, bar_img, tw + gap_px, round(y_bar))
    mpimg.imsave(colorbar_png, combined)
    log(f"Saved colorbar: {colorbar_png}")

    # Overlay: standard "over" alpha compositing onto whichever corner
    # side/vert select (top-right by default), in place -- the combined
    # title+bar image's own transparent margins let the underlying view
    # content show through untouched everywhere else.
    cb_img = mpimg.imread(colorbar_png)
    if side == 'left':
        x0 = side_margin
    else:
        x0 = max(0, vw - side_margin - cb_img.shape[1])
    char_px = 0.6 * label_font_size
    down_shift_px = round(down_shift_chars * char_px)
    if vert == 'top':
        y0 = max(0, vert_margin + down_shift_px)
    else:
        y0 = max(0, vh - vert_margin - cb_img.shape[0] + down_shift_px)
    _alpha_paste(view_img, cb_img, x0, y0)
    mpimg.imsave(output_png, view_img)
    log(f"Overlaid colorbar onto: {output_png}")


def _clip_top_fraction(output_png, frac):
    """Crop off the top `frac` (0-1) of the saved image -- used by
    render_lateral to cut the top 10% of gas headspace (user request,
    2026-08-02). Operates on the raw ParaView render directly (no
    content-based trim happens before this any more, 2026-08-03 -- see
    render_top's/render_lateral's own comments), so `frac` is always a
    fraction of the same fixed height, giving a constant absolute pixel
    crop on every frame."""
    img = mpimg.imread(output_png)
    cut = round(img.shape[0] * frac)
    mpimg.imsave(output_png, img[cut:, :])
    log(f"Clipped top {frac * 100:.0f}%: kept rows [{cut},{img.shape[0] - 1}]")


def _trim_whitespace_bbox(img, pad=8):
    """Crop `img` (RGB or RGBA array) to the bounding box of its non-white
    content on *all four* sides at once. Used for render_transverse's own
    matplotlib output and for domain_schematic.png -- both standalone
    matplotlib figures with wide margins on every side and no fixed
    ParaView-canvas equivalent to rely on instead (unlike
    render_top/render_lateral, which just save the raw fixed-size render
    directly -- see those functions' own comments, 2026-08-03)."""
    non_white = np.any(img[:, :, :3] < 0.98, axis=2)
    rows, cols = np.any(non_white, axis=1), np.any(non_white, axis=0)
    if not rows.any():
        return img
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    h, w = img.shape[:2]
    r0, r1 = max(0, r0 - pad), min(h - 1, r1 + pad)
    c0, c1 = max(0, c0 - pad), min(w - 1, c1 + pad)
    return img[r0:r1 + 1, c0:c1 + 1]


def _trim_sparse_edge(img, side, min_content_px=45, pad=2):
    """Crop away columns at `side` ('left' or 'right') that have fewer than
    min_content_px non-white pixels -- a plain bounding-box trim
    (_trim_whitespace_bbox) still counts a column with just a couple of
    diagonal-edge pixels (an oblique plane's pointed corner) as "content",
    leaving a visually-empty gap at the schematic/cross-section seam even
    though both pieces are already tightly bbox-cropped on their own.
    Translated 1:1 from proto_transverse_3d.py's own trim_sparse_edge (user
    request, 2026-08-03: minimize that gap)."""
    non_white = np.any(img[:, :, :3] < 0.98, axis=2)
    substantial = non_white.sum(axis=0) >= min_content_px
    if not substantial.any():
        return img
    idx = np.where(substantial)[0]
    if side == 'right':
        return img[:, :min(img.shape[1], idx.max() + 1 + pad)]
    return img[:, max(0, idx.min() - pad):]


def _resize_image(img, new_height):
    """Resize `img` (RGB or RGBA array) to exactly `new_height` px, scaling
    width to preserve its own aspect ratio.

    Uses PIL's LANCZOS resampling when PIL is importable -- this repo's own
    lbf3-paraview-mpl image (see Dockerfile.paraview-mpl) has it, pulled in
    as a matplotlib dependency -- for a proper high-quality downscale.
    Falls back to the old "render through a throwaway matplotlib figure"
    trick (imshow + savefig-to-buffer) on the stock kitware/paraview image,
    which has no PIL. That fallback visibly under-resamples a large
    downscale -- e.g. the schematic's ~4800px source rendered down to
    ~400px in the transverse view's composite came out noticeably soft/
    low-quality next to the natively-rendered cross-section content next to
    it (user report, 2026-08-03) -- so prefer the PIL path whenever it's
    available.
    """
    h, w = img.shape[:2]
    new_width = max(1, round(w * new_height / h))
    try:
        from PIL import Image
    except ImportError:
        Image = None
    if Image is not None:
        mode = 'RGBA' if img.shape[2] == 4 else 'RGB'
        pil_img = Image.fromarray((img * 255).astype(np.uint8), mode)
        pil_img = pil_img.resize((new_width, new_height), Image.LANCZOS)
        return np.asarray(pil_img).astype(np.float32) / 255.0
    dpi = 100
    fig = plt.figure(figsize=(new_width / dpi, new_height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy().astype(np.float32) / 255.0
    plt.close(fig)
    return buf


def _draw_cross_section_markers_top(view, y_marker, z_window_min, z_window_max):
    """3 short line markers, one per CROSS_SECTION_X, each spanning
    [z_window_min, z_window_max] at that fixed x, placed at y=y_marker --
    nearer the camera than any real geometry (see render_top's own header),
    so never occluded. Colored via CROSS_SECTION_MARKER_COLORS (dark
    blue/green/red, matching render_transverse's own per-cut colors).

    Callers pass CROSS_SECTION_Z_WINDOW here (not the view's own wider
    Z_VIEW_MIN/MAX crop) -- the lines mark where the transverse view's
    actual cut window sits, not the full top/lateral crop range (user
    request, 2026-08-03: previously spanned the full crop by mistake, a
    leftover from the old distance-behind-laser markers).

    Flat-colored, not scalar-colored: ColorArrayName=['POINTS',''], not
    ColorBy(rep, None) -- the latter crashes when there's no array to
    default to (relevant here since these Line sources carry no data).
    """
    for x_cut, color in zip(CROSS_SECTION_X, CROSS_SECTION_MARKER_COLORS):
        marker = Line(Point1=[x_cut, y_marker, z_window_min], Point2=[x_cut, y_marker, z_window_max])
        disp = Show(marker, view)
        disp.Representation = 'Wireframe'
        disp.ColorArrayName = ['POINTS', '']
        disp.AmbientColor = color
        disp.DiffuseColor = color
        disp.LineWidth = GM_RIM_LINE_WIDTH
        log(f"Cross-section marker: x={x_cut*1e6:+.1f}um -> line at y={y_marker*1e3:.3f}mm, "
            f"z=[{z_window_min*1e3:.3f},{z_window_max*1e3:.3f}]mm")


def _draw_cross_section_frame_lateral(view, x_marker):
    """Single green (x=0 cross-section) rectangle frame spanning
    CROSS_SECTION_Y_WINDOW x CROSS_SECTION_Z_WINDOW, translated to
    x=x_marker -- nearer the camera than any real geometry (see
    render_lateral's own header), so never occluded.

    Only x=0 is drawn, not all 3 cuts: lateral's camera looks down x,
    collapsing every cross-section onto the same (y,z) screen footprint,
    so drawing all 3 frames would just stack 3 identical rectangles
    exactly on top of each other -- the green (middle) one alone conveys
    where the transverse crop sits (user request, 2026-08-03: "show the
    green frame from the lateral view (not the xray)").
    """
    y0, y1 = CROSS_SECTION_Y_WINDOW
    z0, z1 = CROSS_SECTION_Z_WINDOW
    color = CROSS_SECTION_MARKER_COLORS[1]  # middle entry = x=0 = green
    pts = [
        x_marker, y0, z0,
        x_marker, y0, z1,
        x_marker, y1, z1,
        x_marker, y1, z0,
    ]
    frame = PolyLineSource(Points=pts, Closed=1)
    disp = Show(frame, view)
    disp.Representation = 'Wireframe'
    disp.ColorArrayName = ['POINTS', '']
    disp.AmbientColor = color
    disp.DiffuseColor = color
    disp.LineWidth = GM_RIM_LINE_WIDTH
    log(f"Cross-section frame: x=0 -> rect y=[{y0*1e3:.3f},{y1*1e3:.3f}]mm "
        f"z=[{z0*1e3:.3f},{z1*1e3:.3f}]mm at x={x_marker*1e3:.3f}mm")


# ═════════════════════════════════════════════════════════════════════════
# --view=top -- bird's-eye view of the track/spatter/powder-bed layout,
# looking straight down y (build direction). Colored by y (height relative
# to the nominal surface) since that's otherwise invisible from directly
# overhead. Also draws the 3 cross-section marker lines (see
# _draw_cross_section_markers_top): x is the screen-vertical axis here
# (up=-x), so a constant-x cut renders as a horizontal line, drawn at
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
    z_window_min, z_window_max = Z_VIEW_MIN, Z_VIEW_MAX
    log(f"Scan z range=[{laser_table[0][3]*1e3:.3f},{laser_table[-1][3]*1e3:.3f}]mm; "
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
    ycolor.Function = f'(coordsY-{SURFACE_Y})*1e3'  # meters -> mm (was
                                                      # 1e6/um -- user
                                                      # request, 2026-08-02),
                                                      # offset from the
                                                      # surface so 0 = "at
                                                      # the surface",
                                                      # positive = deeper/
                                                      # recessed
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
    ymin_off, ymax_off = (ymin - SURFACE_Y) * 1e3, (ymax - SURFACE_Y) * 1e3
    transition = 0.1  # width of the blue->red transition band (mm, was 100um)
    ctf.RGBPoints = [
        ymin_off,    0.0, 0.0, 1.0,
        -transition, 0.0, 0.0, 1.0,
        transition,  1.0, 0.0, 0.0,
        ymax_off,    1.0, 0.0, 0.0,
    ]

    y_marker = ymin - 0.02 * (ymax - ymin)
    _draw_cross_section_markers_top(view, y_marker, *CROSS_SECTION_Z_WINDOW)

    # Laser center-position marker at the laser's actual current (x, z),
    # same orange as the beam/landing marker in plot_domain_schematic.py
    # for visual consistency across views (user request, 2026-08-02: "show
    # a simple appropriate sign at the laser center position in the top
    # view"). x=0.0: every view in this pipeline already assumes the laser
    # travels along the x=0 centerline (see e.g. plot_domain_schematic.py's
    # TRACK_X), so there's no separate x-position to look up. Placed at the
    # same y_marker plane as the offset lines above -- irrelevant to screen
    # position under this orthographic top-down camera, only used so it
    # isn't accidentally occluded.
    #
    # A hollow orange circle at the laser's true physical spot radius
    # (35um, see CLAUDE.md's "Physical parameters" -- Drude/Fresnel
    # absorptivity, 35 micron radius), replacing a plain filled sphere
    # (user request, 2026-08-02: "draw a circle with the right radius for
    # the laser"). Originally paired with a small white center dot, which
    # was then dropped again -- the circle alone was enough (user
    # follow-up, 2026-08-02).
    # No RegularPolygonSource/Circle/Disk-outline source available in this
    # ParaView build (checked dir(paraview.simple) directly) -- built by
    # hand instead via PolyLineSource, which takes an explicit flat
    # [x0,y0,z0,x1,y1,z1,...] point list and a Closed flag to connect the
    # last point back to the first.
    LASER_SPOT_RADIUS_M = 35e-6
    _n = 48
    _angles = np.linspace(0.0, 2.0 * np.pi, _n, endpoint=False)
    _circle_pts = []
    for _a in _angles:
        _circle_pts.extend([
            LASER_SPOT_RADIUS_M * np.cos(_a), y_marker,
            laser_z + LASER_SPOT_RADIUS_M * np.sin(_a),
        ])
    laser_circle = PolyLineSource(Points=_circle_pts, Closed=1)
    laser_circle_disp = Show(laser_circle, view)
    laser_circle_disp.Representation = 'Wireframe'
    laser_circle_disp.ColorArrayName = ['POINTS', '']
    laser_circle_disp.AmbientColor = [1.0, 0.549, 0.0]
    laser_circle_disp.DiffuseColor = [1.0, 0.549, 0.0]
    laser_circle_disp.LineWidth = 4.0
    # No center dot/x -- the orange circle alone suffices (removed, user
    # request, 2026-08-02).

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

    # No content-based crop here (deliberately, 2026-08-03 -- see this
    # view's own header comment: view.ViewSize above is a function of fixed
    # constants only, so the raw render is already exactly the same pixel
    # size on every timestep. An earlier version cropped this down to
    # whatever pixels happened to be non-blank, which sounds harmless but
    # makes the *saved file's* size track incidental content -- how far the
    # established track/powder bed happens to extend, spatter droplets,
    # where the laser physically is right now -- none of which should be
    # able to change what "the view" is. Traded a small amount of constant
    # blank margin (FRAME_MARGIN already zooms out for this) for dimensions
    # that only ever depend on the fixed crop window, never on the data.
    _overlay_colorbar(ctf, 'y (mm)', output_png, custom_labels=[-0.1, 0.0, 0.1])

    if output_pvsm:
        SaveState(output_pvsm)
        log(f"Saved state: {output_pvsm}")


# ═════════════════════════════════════════════════════════════════════════
# --view=lateral -- through-thickness profile, looking down x. Colored by x
# (otherwise invisible from this orthographic angle) to reveal near/far
# structure. Also draws: (1) a black outline of the mushy (liquid/solid)
# boundary at x=0 (the sample's center depth plane), and (2) a green
# rectangle frame marking where the transverse view's x=0 cross-section
# sits (see _draw_cross_section_frame_lateral -- only x=0 is drawn, since
# this camera collapses all 3 cuts onto the same screen footprint). Both
# are translated to x = xmax + margin -- nearest the camera (camera sits at
# large +x looking toward -x, so larger x is nearer) -- to sit in front of
# the opaque gas/metal surface unoccluded, without changing their (y,z)
# screen position at all (orthographic projection down x never depends on x).
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
    z_window_min, z_window_max = Z_VIEW_MIN, Z_VIEW_MAX
    log(f"Scan z range=[{laser_table[0][3]*1e3:.3f},{laser_table[-1][3]*1e3:.3f}]mm; "
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
    log(f"Mushy-zone (solidus) surface: {ls_poly.GetNumberOfCells()} cells")

    # Liquidus (T=LS_TLIQUIDUS) surface -- same metal_only-then-Contour
    # pattern as the solidus one above, just the other mushy-zone bound.
    ls_liquidus_contour = Contour(Input=metal_only)
    ls_liquidus_contour.ContourBy = ['POINTS', FIELD_LS]
    ls_liquidus_contour.Isosurfaces = [LS_TLIQUIDUS]
    ls_liquidus_contour.UpdatePipeline(time=time_value)
    ls_liquidus_poly = servermanager.Fetch(ls_liquidus_contour)
    log(f"Mushy-zone (liquidus) surface: {ls_liquidus_poly.GetNumberOfCells()} cells")

    def _outline_at_x0(contour):
        """Exact plane intersection at x=0, cropped to the view window --
        not a rendered clip, so there's no depth-buffer ambiguity to
        arbitrate for this 1D curve (same technique as
        render_transverse's slice_at_cut()). Shared by the solidus and
        liquidus outlines below (identical pipeline, different contour
        field value feeding in).

        Bug fix (solidus case originally): `contour` is built off the
        *whole* domain (unclipped, unlike gm_contour -> feature below) --
        the x=0 slice through it can carry mushy-zone content at any z/y
        across the entire mesh, not just inside
        [z_window_min,z_window_max] x [Y_DEPTH_MIN,Y_DEPTH_MAX]. Camera
        framing alone doesn't stop this from showing up:
        FRAME_MARGIN's zoom-out margin means the actual rendered z-range
        is wider than [z_window_min,z_window_max] (by the same factor),
        so real mushy-outline content just past the intended window was
        visible in the saved image, undermining the whole point of
        fixing Z_VIEW_MIN/MAX. Explicit box-clip, matching feature's own
        bounds.
        """
        s = Slice(Input=contour)
        s.SliceType = 'Plane'
        s.SliceType.Origin = [0.0, (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0, (z_window_min + z_window_max) / 2.0]
        s.SliceType.Normal = [1.0, 0.0, 0.0]
        s.UpdatePipeline(time=time_value)
        c = Clip(Input=s)
        c.ClipType = 'Box'
        c.ClipType.Position = [-1e-6, Y_DEPTH_MIN, z_window_min]
        c.ClipType.Length = [2e-6, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
        c.Invert = 1
        c.UpdatePipeline(time=time_value)
        return c

    ls_slice_clip = _outline_at_x0(ls_contour)
    ls_slice_poly = servermanager.Fetch(ls_slice_clip)
    log(f"Solidus outline at x=0 (clipped to view window): {ls_slice_poly.GetNumberOfCells()} cells "
        f"(may be empty if no melt currently straddles the sample's center depth)")

    ls_liquidus_slice_clip = _outline_at_x0(ls_liquidus_contour)
    ls_liquidus_slice_poly = servermanager.Fetch(ls_liquidus_slice_clip)
    log(f"Liquidus outline at x=0 (clipped to view window): {ls_liquidus_slice_poly.GetNumberOfCells()} cells "
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
    xcolor.Function = 'coordsX*1e3'  # meters -> mm (was 1e6/um -- user request, 2026-08-02)
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
    xmin_mm, xmax_mm = xmin * 1e3, xmax * 1e3
    transition = 0.05  # width of the blue->red transition band (mm, was 50um)
    ctf.RGBPoints = [
        xmin_mm,     0.0, 0.0, 1.0,
        -transition, 0.0, 0.0, 1.0,
        transition,  1.0, 0.0, 0.0,
        xmax_mm,     1.0, 0.0, 0.0,
    ]

    # Shared "front of camera" x position for both the mushy outline and
    # the transverse cut markers below -- see this function's header.
    x_marker = xmax + 0.02 * (xmax - xmin)

    ls_outline = Transform(Input=ls_slice_clip)
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

    # Liquidus outline -- same x=0/x_marker translate as the solidus one
    # above, fainter gray so the two read as distinct at a glance (liquidus
    # nested inside solidus, since anything past liquidus is already fully
    # solidus too).
    ls_liquidus_outline = Transform(Input=ls_liquidus_slice_clip)
    ls_liquidus_outline.Transform = 'Transform'
    ls_liquidus_outline.Transform.Translate = [x_marker, 0.0, 0.0]
    ls_liquidus_outline.UpdatePipeline(time=time_value)
    ls_liquidus_outline_disp = Show(ls_liquidus_outline, view)
    ls_liquidus_outline_disp.Representation = 'Wireframe'
    ls_liquidus_outline_disp.ColorArrayName = ['POINTS', '']
    ls_liquidus_outline_disp.AmbientColor = LS_LIQUIDUS_COLOR
    ls_liquidus_outline_disp.DiffuseColor = LS_LIQUIDUS_COLOR
    ls_liquidus_outline_disp.LineWidth = LS_LINE_WIDTH

    _draw_cross_section_frame_lateral(view, x_marker)

    # Camera: looking down +x, up = -y so atmosphere is "up" in the frame.
    view.CameraParallelProjection = 1
    view.CameraViewUp = [0, -1, 0]
    view.CameraFocalPoint = [x_center, y_center, z_center]
    view.CameraPosition = [x_center + 2.0 * (xmax - xmin), y_center, z_center]
    view.CameraParallelScale = (Y_DEPTH_MAX - Y_DEPTH_MIN) / 2.0 * FRAME_MARGIN
    Render(view)

    SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
    log(f"Saved: {output_png}")

    # No content-based crop here (see render_top's own comment on the same
    # change, 2026-08-03) -- view.ViewSize is a function of fixed constants
    # only, so skipping it means the saved file's dimensions can never
    # depend on melt-pool depth, spatter, or where the laser physically is.
    # _clip_top_fraction below is unaffected by this -- it's already a
    # fixed *fraction* of whatever height gets passed in, and that height
    # is now itself constant, so the absolute pixel crop it applies is
    # constant too (previously it wasn't, since it ran on the
    # already-content-cropped, variable-height image).
    _clip_top_fraction(output_png, 0.10)
    _overlay_colorbar(ctf, 'x (mm)', output_png, custom_labels=[-0.1, 0.0, 0.1])

    if output_pvsm:
        SaveState(output_pvsm)
        log(f"Saved state: {output_pvsm}")


# ═════════════════════════════════════════════════════════════════════════
# --view=transverse -- 3 cross-section cuts (X_CROSS_SECTIONS) shown as a
# single oblique 3D orthographic scene, not the old flat same-(y,z)-
# footprint overlay (replaced 2026-08-03, user-driven redesign): each cut is
# its own plane in space, spread apart by 3 small sign-based +-offsets (x,
# scan-direction z, and depth y -- see X_STAGGER_STEP_UM/Z_STAGGER_STEP_MM/
# Y_STAGGER_STEP_UM below) instead of the one big multiplicative x
# exaggeration tried earlier. Each plane carries its own gas/metal boundary
# curve (color-coded, same red/green/blue as the old view) and a semi-
# transparent fill for the metal (non-gas) side, itself split at the
# solidus temperature into solid (the cut's own base color) and liquid (a
# darker shade of it, closer to the gas/metal curve's own full-saturation
# color) -- plus a dotted white liquidus outline on top (all 3, user
# request, 2026-08-03). The gas side of the window is left empty. A dashed
# "ghost" frame at each cut's TRUE (unoffset) position, plus a thin leader
# line to it, shows where it actually sits so the offsets don't silently
# mislead (no numeric callout in the final frames -- agreed with the user
# to keep that honesty check purely visual; skipped for the middle/x=0 cut,
# which carries no offset at all). All of this was tuned
# interactively against proto_transverse_3d.py's synthetic sin curves (much
# faster than this real-data Docker pipeline) and then translated back here
# 1:1.
#
# Rendered via a hand-rolled painter's algorithm on a plain 2D matplotlib
# axes (mplot3d's own 3D rendering turned out to mismatch cross-plane draw
# order too often -- see the function's own header comment further down for
# why), not ParaView's own Show()/Render() -- ParaView is used purely to
# extract the sliced geometry (same "ParaView extracts, matplotlib draws"
# split render_xray already uses for its own, unrelated ray-tracing
# technique).
# Prototyped first with synthetic sin curves in proto_transverse_3d.py
# before being wired to this real data.
# ═════════════════════════════════════════════════════════════════════════
def render_transverse(foam_file, time_value, output_png, output_pvsm):
    # Cross-section x-positions (m) -- unchanged from the old view.
    X_CROSS_SECTIONS = [-25e-6, 0.0, 25e-6]
    # 3 explicit shades per hue (user-specified, 2026-08-03), replacing the
    # earlier _darken_hex()-computed shades -- negative cross-section is
    # blue, positive is red (flipped, user request, 2026-08-03; was
    # red/green/blue), translated from proto_transverse_3d.py. LIGHT is the
    # solid-metal fill color (was CROSS_SECTION_COLORS), MEDIUM the liquid
    # fill, DARK the solid/ghost frame borders, leader lines, and gas/metal
    # outline.
    LIGHT_COLORS = ['#2a78d6', '#66bb6a', '#e34948']   # light blue/green/red
    MEDIUM_COLORS = ['#3f51b5', '#1f8b4d', '#d32f2f']  # medium blue/green/red
    DARK_COLORS = ['#0d47a1', '#256428', '#b71c1c']    # dark blue/green/red
                                     # -- green darkened further in both
                                     # (was #27ae60/#2e7d32), user request,
                                     # 2026-08-03
    CROSS_SECTION_COLORS = LIGHT_COLORS  # kept as an alias -- used below
    # z window (m), fixed in absolute z (not laser-relative) so cross-
    # sections stay comparable across an entire batch of timesteps.
    Z_WINDOW_MIN, Z_WINDOW_MAX = 1.5e-3, 2.0e-3
    # y window (m), shared across all cross-sections.
    Y_MIN, Y_MAX = 0.225e-3, 0.35e-3

    # Camera angles -- tuned together with the user against
    # proto_transverse_3d.py; AZIMUTH_DEG is named (not inlined into
    # view_init) so it stays the one thing to point at when discussing/
    # adjusting the view's left-right rotation. No horizontal (x)
    # exaggeration anymore -- tried up to 50x (2026-08-03) combined with
    # true x:z scale, which made the scene far too wide and shrank the
    # planes to illegible slivers; replaced by the vertical y-stagger
    # below, which the user asked for specifically so the 3 planes read as
    # directly stacked ("right on top of each other") rather than spread
    # out sideways.
    # Matches plot_domain_schematic.py's own view_init(elev=20, azim=-40)
    # exactly (was elev=18/azim=-60) -- user request, 2026-08-03: this view
    # should read as a zoomed-in version of the schematic's own small
    # cross-section marker frames, which only works if both use the same
    # camera angle.
    ELEVATION_DEG = 20
    AZIMUTH_DEG = 40

    SCHEMATIC_GAP_PX = 0                   # white gap between the schematic
                                            # and the overlay
    SCHEMATIC_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     'plot_domain_schematic.py')  # regenerated
                                     # per-frame below (laser position/time
                                     # label vary), not a static pre-baked
                                     # image

    reader = OpenFOAMReader(FileName=foam_file)
    reader.CellArrays = [FIELD_GM, FIELD_LS]  # FIELD_LS back for the
                                               # solidus/liquidus split
                                               # below (was gas/metal only)
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

    z_window_min, z_window_max = Z_WINDOW_MIN, Z_WINDOW_MAX
    log(f"Fixed z window: z=[{z_window_min*1e3:.3f},{z_window_max*1e3:.3f}]mm; "
        f"y=[{Y_MIN*1e3:.3f},{Y_MAX*1e3:.3f}]mm")

    gm_contour = Contour(Input=merged)
    gm_contour.ContourBy = ['POINTS', FIELD_GM]
    gm_contour.Isosurfaces = [ISO_THRESHOLD]
    gm_contour.UpdatePipeline(time=time_value)
    gm_poly = servermanager.Fetch(gm_contour)
    log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

    # Liquidus contour (dotted white overlay, user request, 2026-08-03),
    # restricted to inside the metal first -- same Clip-then-Contour
    # pattern render_lateral uses (see its own comment on why Clip, not
    # Threshold). Solidus itself is never drawn as a curve -- it's only
    # used below (per x-cut) to split the metal fill into solid/liquid.
    metal_only = Clip(Input=merged)
    metal_only.ClipType = None
    metal_only.Scalars = ['POINTS', FIELD_GM]
    metal_only.Value = ISO_THRESHOLD
    metal_only.Invert = 0
    metal_only.UpdatePipeline(time=time_value)

    ls_liquidus_contour = Contour(Input=metal_only)
    ls_liquidus_contour.ContourBy = ['POINTS', FIELD_LS]
    ls_liquidus_contour.Isosurfaces = [LS_TLIQUIDUS]
    ls_liquidus_contour.UpdatePipeline(time=time_value)
    ls_liquidus_poly = servermanager.Fetch(ls_liquidus_contour)
    log(f"Liquidus surface: {ls_liquidus_poly.GetNumberOfCells()} cells")

    y_center = (Y_MIN + Y_MAX) / 2.0
    z_center = (z_window_min + z_window_max) / 2.0

    def box_clip(input_poly, x_window_min, x_window_max):
        """y/z fixed crop window (Y_MIN/MAX, z_window_min/max), x in
        [x_window_min, x_window_max]."""
        c = Clip(Input=input_poly)
        c.ClipType = 'Box'
        c.ClipType.Position = [x_window_min, Y_MIN, z_window_min]
        c.ClipType.Length = [x_window_max - x_window_min, Y_MAX - Y_MIN, z_window_max - z_window_min]
        c.Invert = 1
        c.UpdatePipeline(time=time_value)
        return c

    def slice_at_x(input_poly, x_cut):
        """Exact plane intersection at x=x_cut, then crop to the y/z window
        -- a Slice computes the exact 1D intersection curve mathematically
        (same technique as render_lateral's own x=0 outline slices)."""
        s = Slice(Input=input_poly)
        s.SliceType = 'Plane'
        s.SliceType.Origin = [x_cut, y_center, z_center]
        s.SliceType.Normal = [1.0, 0.0, 0.0]
        s.UpdatePipeline(time=time_value)
        return box_clip(s, x_cut - 1e-9, x_cut + 1e-9)

    # Cell-type-filtered GetCellPoints() traversal, not GetLines()/GetPolys()
    # (the technique _load_laser_rays uses) -- these two helpers' inputs are
    # Clip outputs, always vtkUnstructuredGrid (never vtkPolyData, which is
    # the only class GetLines()/GetPolys() exist on), so traversal has to go
    # through the dataset-agnostic GetCellType()/GetCellPoints() API instead.
    def extract_polylines(vtk_ds):
        """List of (N,3) point arrays, one per line-type cell."""
        if vtk_ds.GetNumberOfPoints() == 0:
            return []
        pts = vtk_to_numpy(vtk_ds.GetPoints().GetData())
        id_list = vtk.vtkIdList()
        segments = []
        for i in range(vtk_ds.GetNumberOfCells()):
            if vtk_ds.GetCellType(i) not in (vtk.VTK_LINE, vtk.VTK_POLY_LINE):
                continue
            vtk_ds.GetCellPoints(i, id_list)
            n = id_list.GetNumberOfIds()
            segments.append(np.array([pts[id_list.GetId(k)] for k in range(n)]))
        return segments

    def extract_polygons(vtk_ds):
        """List of (N,3) point arrays, one per polygon-type cell -- vertex
        loops ready for projection/fill."""
        if vtk_ds.GetNumberOfPoints() == 0:
            return []
        pts = vtk_to_numpy(vtk_ds.GetPoints().GetData())
        id_list = vtk.vtkIdList()
        cells = []
        poly_types = (vtk.VTK_TRIANGLE, vtk.VTK_QUAD, vtk.VTK_POLYGON)
        for i in range(vtk_ds.GetNumberOfCells()):
            if vtk_ds.GetCellType(i) not in poly_types:
                continue
            vtk_ds.GetCellPoints(i, id_list)
            n = id_list.GetNumberOfIds()
            cells.append(np.array([pts[id_list.GetId(k)] for k in range(n)]))
        return cells

    tmp_dir = tempfile.mkdtemp(prefix="transverse_")

    try:
        # This view no longer renders through mplot3d's own draw path at
        # all. mplot3d's automatic depth-sort (default) and its manual
        # mode (computed_zorder=False + explicit per-artist zorder, the
        # same fix plot_domain_schematic.py uses for its own "line behind
        # glass" issue) *both* turned out to mix up the draw order between
        # different cross-section planes' overlapping fills/lines, even
        # with widely-separated zorder values -- a real mplot3d limitation
        # with many overlapping 3D artists, not a sign/logic bug on our end
        # (the near/far assignment itself was verified correct via
        # proj3d.proj_transform's own depth output). User report,
        # 2026-08-03: "the green outline is not at appropriate order it is
        # shown on top ... the +0.25 cross-section".
        #
        # Fix: a manual painter's algorithm. `ax` below is never drawn/
        # saved -- it exists purely so get_proj() can hand back the exact
        # projection matrix (M) matplotlib would use for this view (same
        # elev/azim/box_aspect/limits as before). Every primitive (each
        # small fill polygon, each curve segment, each frame/leader line,
        # each label) is projected to 2D via that M *by hand*
        # (proj3d.proj_transform -- the same function
        # plot_domain_schematic.py's own callout-placement code already
        # uses) right when it's built, carrying its own camera-depth
        # (smaller = nearer, verified empirically). Once every plane's
        # primitives are collected, they're sorted by that depth and drawn
        # on a plain 2D axes in far-to-near order -- a 2D Axes reliably
        # draws artists in call order, so this sidesteps mplot3d's
        # unreliable cross-artist depth sorting entirely.
        #
        # Figsize aspect (30:22) still matches plot_domain_schematic.py's
        # own figsize exactly -- the pixel-box aspect ratio that matters
        # for the schematic's green frame and this plot's own plane edges
        # to project at matching on-screen slopes is baked into M itself.
        fig3d = plt.figure(figsize=(10, 10 * 22 / 30))
        ax = fig3d.add_axes([0, 0, 1, 1], projection='3d')
        ax.set_proj_type('ortho')  # no perspective

        def rect_lines(x_um, z_off_mm=0.0, y_off_um=0.0):
            """The y=[Y_MIN,Y_MAX] x z=[z_window_min,z_window_max] window-
            frame rectangle at a given (already-displayed-units) x, shifted
            by z_off_mm/y_off_um, as (xs, zs_mm, ys_um) plot-ready lists --
            matches proto_transverse_3d.py's own rect_lines exactly."""
            rect_z_mm = [z_window_min * 1e3 + z_off_mm, z_window_max * 1e3 + z_off_mm,
                         z_window_max * 1e3 + z_off_mm, z_window_min * 1e3 + z_off_mm,
                         z_window_min * 1e3 + z_off_mm]
            rect_y_um = [Y_MIN * 1e6 + y_off_um, Y_MIN * 1e6 + y_off_um,
                         Y_MAX * 1e6 + y_off_um, Y_MAX * 1e6 + y_off_um, Y_MIN * 1e6 + y_off_um]
            return [x_um] * 5, rect_z_mm, rect_y_um

        # 3 artificial offsets, translated 1:1 from proto_transverse_3d.py
        # (tuned there interactively -- much faster than the real-data
        # Docker pipeline -- then ported back here, user request,
        # 2026-08-03: "let's translate back to our actual transverse
        # cross-sections"). All 3 are sign-based +-steps off the middle
        # (x=0) cross-section, which stays fully at its true position in
        # x/y/z alike:
        #  - x: +-137.5um, left/right
        #  - z (scan direction): +-135um, back/forward
        #  - y (depth): +-6.75um, up/down
        # Each cut's ghost frame (true position) and leader line are
        # skipped entirely when all 3 offsets are zero (the middle cut) --
        # it would just exactly coincide with the solid frame.
        X_STAGGER_STEP_UM = 137.5
        Z_STAGGER_STEP_MM = 0.135
        Y_STAGGER_STEP_UM = 6.75
        x_offsets_um = [X_STAGGER_STEP_UM if x > 0 else (-X_STAGGER_STEP_UM if x < 0 else 0.0)
                        for x in X_CROSS_SECTIONS]
        z_offsets_mm = [Z_STAGGER_STEP_MM if x > 0 else (-Z_STAGGER_STEP_MM if x < 0 else 0.0)
                        for x in X_CROSS_SECTIONS]
        y_offsets_um = [Y_STAGGER_STEP_UM if x > 0 else (-Y_STAGGER_STEP_UM if x < 0 else 0.0)
                        for x in X_CROSS_SECTIONS]

        def apply_box_aspect(x0, x1, y0, y1, z0, z1, aspect):
            """Set the 3 axes' view limits to render in the given relative
            proportions -- same pre-3.3 mplot3d fallback as
            plot_domain_schematic.py's own apply_box_aspect (that Docker
            image's matplotlib 3.1.1 has no set_box_aspect). Called with
            *ascending* z0<z1 -- ax.invert_zaxis() below handles the
            deeper-points-down flip separately, same order the schematic
            uses, since inverting first would make z1-z0 negative and break
            the span math."""
            if hasattr(ax, 'set_box_aspect'):
                ax.set_xlim3d(x0, x1)
                ax.set_ylim3d(y0, y1)
                ax.set_zlim3d(z0, z1)
                ax.set_box_aspect(aspect)
                return
            true_ranges = (x1 - x0, y1 - y0, z1 - z0)
            centers = ((x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2)
            max_r = max(aspect)
            spans = [t * max_r / r for t, r in zip(true_ranges, aspect)]
            (sx0, sx1), (sy0, sy1), (sz0, sz1) = (
                (c - s / 2, c + s / 2) for c, s in zip(centers, spans))
            ax.set_xlim3d(sx0, sx1)
            ax.set_ylim3d(sy0, sy1)
            ax.set_zlim3d(sz0, sz1)

        # Axis limits cover both the displayed (offset) solid frames and the
        # true (unoffset) ghost frames, on all 3 axes now.
        xs_disp_um = [x * 1e6 + off for x, off in zip(X_CROSS_SECTIONS, x_offsets_um)]
        x_lo, x_hi = min(xs_disp_um + [x * 1e6 for x in X_CROSS_SECTIONS]) - 20, \
            max(xs_disp_um + [x * 1e6 for x in X_CROSS_SECTIONS]) + 20
        z_candidates_mm = (
            [z_window_min * 1e3 + off for off in z_offsets_mm]
            + [z_window_max * 1e3 + off for off in z_offsets_mm]
            + [z_window_min * 1e3, z_window_max * 1e3]
        )
        z_lo_mm, z_hi_mm = min(z_candidates_mm), max(z_candidates_mm)
        y_candidates_um = (
            [Y_MIN * 1e6 + off for off in y_offsets_um]
            + [Y_MAX * 1e6 + off for off in y_offsets_um]
            + [Y_MIN * 1e6, Y_MAX * 1e6]
        )
        y_lo_um, y_hi_um = min(y_candidates_um), max(y_candidates_um)
        # Same real-world scale for x and z (user request, 2026-08-03: the
        # rendered cross-section rectangles looked square instead of long,
        # "not intuitive ... let's have z and x direction to have the same
        # size") -- aspect components are the axes' *true* ranges, all
        # converted to one consistent unit (um) even though the axis limits
        # themselves stay in the plot's own mixed display units (um for
        # x/y, mm for z). This also restores the rectangles' true (long,
        # z-span 4x the y-span) shape, matching plot_domain_schematic.py's
        # own inset.
        apply_box_aspect(x_lo, x_hi, z_lo_mm, z_hi_mm, y_lo_um, y_hi_um,
                          (x_hi - x_lo, (z_hi_mm - z_lo_mm) * 1e3, y_hi_um - y_lo_um))
        ax.invert_zaxis()  # deeper (larger real-y) points down visually
        ax.view_init(elev=ELEVATION_DEG, azim=-AZIMUTH_DEG)
        fig3d.canvas.draw()  # settle the projection before reading it back
        M = ax.get_proj()
        plt.close(fig3d)  # never rendered/saved -- only used to compute M

        # Every primitive is (kind, x2d, y2d, depth, draw_kwargs); collected
        # here, drawn far-to-near (largest depth first) further below.
        primitives = []

        def add_line(xs, ys, zs, **kwargs):
            x2, y2, z2 = proj3d.proj_transform(np.asarray(xs, dtype=float),
                                                np.asarray(ys, dtype=float),
                                                np.asarray(zs, dtype=float), M)
            primitives.append(('line', x2, y2, float(np.mean(z2)), kwargs))

        def add_fill(xs, ys, zs, **kwargs):
            x2, y2, z2 = proj3d.proj_transform(np.asarray(xs, dtype=float),
                                                np.asarray(ys, dtype=float),
                                                np.asarray(zs, dtype=float), M)
            primitives.append(('fill', x2, y2, float(np.mean(z2)), kwargs))

        def add_text(x, y, z, s, dx_pt, dy_pt, **kwargs):
            x2, y2, z2 = proj3d.proj_transform(np.array([x]), np.array([y]), np.array([z]), M)
            kwargs.update(s=s, dx_pt=dx_pt, dy_pt=dy_pt)
            primitives.append(('text', x2, y2, float(z2[0]), kwargs))

        for x_cut, color, medium_color, dark_color, x_off_um, z_off_mm, y_off_um in zip(
                X_CROSS_SECTIONS, LIGHT_COLORS, MEDIUM_COLORS, DARK_COLORS,
                x_offsets_um, z_offsets_mm, y_offsets_um):
            x_true_um = x_cut * 1e6
            x_disp_um = x_true_um + x_off_um
            log(f"x={x_true_um:.1f}um (color {color}) -> displayed at {x_disp_um:.1f}um, "
                f"z-offset={z_off_mm*1e3:+.1f}um, y-offset={y_off_um:+.1f}um")

            gm_rim = slice_at_x(gm_contour, x_cut)
            gm_rim_poly = servermanager.Fetch(gm_rim)
            curve_segments = extract_polylines(gm_rim_poly)
            log(f"  Gas/metal outline: {len(curve_segments)} segments "
                f"({gm_rim_poly.GetNumberOfCells()} cells)")
            for seg in curve_segments:
                add_line([x_disp_um] * len(seg), seg[:, 2] * 1e3 + z_off_mm, seg[:, 1] * 1e6 + y_off_um,
                         color=dark_color, linewidth=1.0)

            # Liquidus outline data -- extracted and logged, but not
            # currently drawn (dotted white overlay disabled for now, see
            # below -- user request, 2026-08-03).
            liquidus_rim = slice_at_x(ls_liquidus_contour, x_cut)
            liquidus_rim_poly = servermanager.Fetch(liquidus_rim)
            liquidus_segments = extract_polylines(liquidus_rim_poly)
            log(f"  Liquidus outline: {len(liquidus_segments)} segments "
                f"(may be empty if no liquid straddles this x/window)")

            # Metal-region fill, split into solid and liquid (user request,
            # 2026-08-03) -- slice the full field (not just the contour),
            # keep the metal side (Invert=0, same convention as before),
            # then further split on FIELD_LS at the solidus temperature:
            # Invert=0 keeps T>=solidus (liquid, including the mushy zone),
            # Invert=1 keeps T<solidus (solid) -- same "Invert=0 keeps
            # values above the threshold" convention already established by
            # the gas/metal clip above. Each half cropped to the y/z window
            # and filled separately; left empty if this cut has no
            # solid/liquid metal in the window at all.
            field_slice = Slice(Input=merged)
            field_slice.SliceType = 'Plane'
            field_slice.SliceType.Origin = [x_cut, y_center, z_center]
            field_slice.SliceType.Normal = [1.0, 0.0, 0.0]
            field_slice.UpdatePipeline(time=time_value)

            metal_clip = Clip(Input=field_slice)
            metal_clip.ClipType = None
            metal_clip.Scalars = ['POINTS', FIELD_GM]
            metal_clip.Value = ISO_THRESHOLD
            metal_clip.Invert = 0
            metal_clip.UpdatePipeline(time=time_value)

            def clip_temp(input_poly, value, invert):
                """One scalar Clip on FIELD_LS -- Invert=0 keeps T>=value,
                Invert=1 keeps T<value (same convention as the gas/metal
                clip above)."""
                c = Clip(Input=input_poly)
                c.ClipType = None
                c.Scalars = ['POINTS', FIELD_LS]
                c.Value = value
                c.Invert = invert
                c.UpdatePipeline(time=time_value)
                return c

            # 3-way split -- solid (T<solidus, light shade), mushy
            # (solidus<=T<liquidus, dark shade -- same shade as the
            # borders/outline, user request, 2026-08-03: "use dark color to
            # color the part of the liquid region that is between liquidus
            # and solidus"), and fully liquid (T>=liquidus, medium shade).
            # mushy is 2 chained clips (T>=solidus, then T<liquidus on that
            # result); solid/liquid are each a single clip, same convention
            # as the gas/metal clip above.
            FILL_ALPHA_SOLID = 1.0  # was 0.5, 0.8, 1.0, 0.30, 0.35 before
            FILL_ALPHA_MUSHY = 1.0
            FILL_ALPHA_LIQUID = 1.0  # -- user request, 2026-08-03: full opacity
            for phase_name, phase_color, phase_alpha, phase_result in (
                ('solid', color, FILL_ALPHA_SOLID, clip_temp(metal_clip, LS_TSOLIDUS, 1)),
                ('mushy', dark_color, FILL_ALPHA_MUSHY,
                 clip_temp(clip_temp(metal_clip, LS_TSOLIDUS, 0), LS_TLIQUIDUS, 1)),
                ('liquid', medium_color, FILL_ALPHA_LIQUID, clip_temp(metal_clip, LS_TLIQUIDUS, 0)),
            ):
                phase_boxed = box_clip(phase_result, x_cut - 1e-9, x_cut + 1e-9)
                phase_poly = servermanager.Fetch(phase_boxed)
                fill_cells = extract_polygons(phase_poly)
                log(f"  {phase_name.capitalize()} fill: {len(fill_cells)} polygons "
                    f"(may be empty if no {phase_name} metal straddles this x/window)")
                # One add_fill() per mesh cell (not one big collection for
                # the whole phase) -- each cell needs its own depth for the
                # cross-plane painter's algorithm above to work at all.
                for cell in fill_cells:
                    add_fill([x_disp_um] * len(cell), [p[2] * 1e3 + z_off_mm for p in cell],
                             [p[1] * 1e6 + y_off_um for p in cell],
                             facecolor=phase_color, edgecolor=phase_color,
                             alpha=phase_alpha, linewidth=0.3)

            # Dotted white liquidus overlay disabled for now (user request,
            # 2026-08-03) -- liquidus_segments is still fetched/logged
            # above so the data's there to turn this back on later.

            # Solid cross-section window frame at the displayed x/z/y --
            # dark shade (user request, 2026-08-03), like the ghost frame/
            # leader line/gas outline below.
            rx, rz, ry = rect_lines(x_disp_um, z_off_mm, y_off_um)
            add_line(rx, rz, ry, color=dark_color, linewidth=1.0, alpha=0.6)

            # Dashed ghost frame at the true (unoffset) x/z/y, plus a thin
            # leader line (top-far corner) connecting it to the displayed
            # position -- skipped entirely when there's no offset at all
            # (the middle/x=0 cross-section), since the ghost frame would
            # just exactly coincide with the solid one.
            if x_off_um != 0 or z_off_mm != 0 or y_off_um != 0:
                gx, gz, gy = rect_lines(x_true_um)
                add_line(gx, gz, gy, color=dark_color, linewidth=1.2, linestyle=(0, (3, 2)), alpha=0.9)
                add_line([x_true_um, x_disp_um], [z_window_max * 1e3, z_window_max * 1e3 + z_off_mm],
                         [Y_MIN * 1e6, Y_MIN * 1e6 + y_off_um],
                         color=dark_color, linewidth=0.8, linestyle=(0, (1, 2)), alpha=0.7)

                # Legend-style label at the ghost frame's bottom-left corner
                # (low z, high/deep y), same color as the frame, just
                # outside it, nudged in screen-space (points) -- translated
                # 1:1 from proto_transverse_3d.py.
                label_fontsize = 5.5 * 1.2 * 1.25 * 1.3 * 1.3  # bumped again, user request, 2026-08-03 (x2)
                char_pt = 0.6 * label_fontsize
                # Proper Unicode minus (U+2212), not the plain ASCII hyphen.
                val = x_true_um / 100
                sign = '+' if val > 0 else ('−' if val < 0 else '')
                label = f"x={sign}{abs(val):.2f} μm"  # user request, 2026-08-03 -- verified
                                                       # this Docker image's default font
                                                       # (DejaVu Sans) does render μ correctly
                # dx/dy in character units (positive dx = right, positive dy = up).
                dx_chars = -10 if x_true_um < 0 else 9  # positive nudged 2 more chars left, user request, 2026-08-03
                dy_chars = 5 if x_true_um < 0 else 3  # negative nudged 1 char down, positive nudged 2 more chars down, user request, 2026-08-03
                add_text(x_true_um, z_window_min * 1e3 - 0.03, Y_MAX * 1e6 + 20, label,
                         dx_chars * char_pt, dy_chars * char_pt,
                         color=dark_color, fontsize=label_fontsize, ha='left', va='top')

        # Draw everything on a plain 2D axes, farthest-to-nearest -- a 2D
        # Axes reliably draws artists in call order, so this is the actual
        # painter's algorithm (see this function's header comment above for
        # why that's done by hand instead of trusting mplot3d's own
        # sorting). set_aspect('equal') keeps the proportions M already
        # computed (box_aspect included) instead of independently
        # re-stretching x vs y to fill the axes box, which would undo the
        # whole point.
        fig = plt.figure(figsize=(10, 10 * 22 / 30))
        ax2d = fig.add_axes([0.02, 0.02, 0.96, 0.96])
        ax2d.set_aspect('equal')
        ax2d.axis('off')
        for kind, x2, y2, depth, kwargs in sorted(primitives, key=lambda p: p[3], reverse=True):
            if kind == 'fill':
                ax2d.fill(x2, y2, **kwargs)
            elif kind == 'line':
                ax2d.plot(x2, y2, **kwargs)
            elif kind == 'text':
                text_kwargs = {k: v for k, v in kwargs.items() if k not in ('s', 'dx_pt', 'dy_pt')}
                txt = ax2d.text(x2[0], y2[0], kwargs['s'], **text_kwargs)
                offset = offset_copy(txt.get_transform(), fig=fig,
                                      x=kwargs['dx_pt'], y=kwargs['dy_pt'], units='points')
                txt.set_transform(offset)
        ax2d.relim()
        ax2d.autoscale_view()

        # dpi bumped from 150 -- at the same dpi as proto_transverse_3d.py,
        # this still came out lower-resolution (1149x240 vs proto's
        # 1425x382, user report, 2026-08-03): the real gas/metal curve
        # doesn't span the full z window edge-to-edge the way the synthetic
        # sin curve does by construction, so its content bounding box (what
        # gets kept after the crop below) occupies a smaller fraction of the
        # same-sized rendered canvas. Bumped to compensate.
        fig.savefig(output_png, dpi=220, bbox_inches='tight')
        plt.close(fig)
        # bbox_inches='tight' alone leaves a lot of dead margin here --
        # explicit content-based crop (same helper plot_domain_schematic.py
        # uses for its own wide-margin figure) closes it.
        mpimg.imsave(output_png, _trim_whitespace_bbox(mpimg.imread(output_png), pad=15))
        log(f"Saved: {output_png}")

        if output_pvsm:
            log("Skipping .pvsm save for --view=transverse -- matplotlib-rendered, "
                "no ParaView RenderView/state involved (like --view=xray)")

        # Prepend the domain schematic to the left of the overlay.
        # Regenerated fresh for *this* frame (not a static pre-baked image)
        # -- its single laser marker/time label must match the actual
        # timestep being rendered here, so plot_domain_schematic.py is
        # invoked as a subprocess with this frame's own laser_z and a time
        # label, output to tmp_dir (cleaned up automatically below).
        # Cropped to its own content bounding box (its source figure has
        # wide margins on every side) and resized to the row's own height,
        # preserving aspect ratio.
        if os.path.exists(SCHEMATIC_SCRIPT):
            schematic_png = os.path.join(tmp_dir, 'schematic.png')
            # Plain ASCII number, not "... μs" -- passed as a subprocess
            # argv element below, and this Docker image's pvpython (Python
            # 3.6) encodes argv using the *parent* process's own locale
            # (fixed at interpreter startup, unaffected by an env= passed to
            # subprocess.run), which is ASCII/POSIX here, so a μ character
            # in an argv string raises UnicodeEncodeError before the child
            # process even starts. plot_domain_schematic.py appends its own
            # "μs" suffix internally instead (literal in its own source, no
            # argv involved) -- user request, 2026-08-03.
            time_label = f"{time_value * 1e6:.0f}"
            try:
                # sys.executable resolves to an internal VTK launcher helper
                # here (not directly invocable), not the real pvpython
                # binary -- use the same absolute path every shell script in
                # this repo already hardcodes for this Docker image instead.
                subprocess.run(
                    ['/opt/paraview/bin/pvpython', SCHEMATIC_SCRIPT,
                     f"{laser_z * 1e6:.3f}", time_label, schematic_png],
                    check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    universal_newlines=True,  # 'text=True' equivalent -- this
                                               # Docker image's pvpython is
                                               # Python 3.6, predates 'text='
                )
            except subprocess.CalledProcessError as e:
                log(f"plot_domain_schematic.py failed (skipping schematic): {e.stderr}")
                schematic_png = None
            if schematic_png and os.path.exists(schematic_png):
                row_img = _trim_sparse_edge(mpimg.imread(output_png), 'left')
                schematic = _trim_sparse_edge(_trim_whitespace_bbox(mpimg.imread(schematic_png)), 'right')
                # Schematic occupies exactly 25% of the final combined
                # width, the cross-section render 75% (user request,
                # 2026-08-03) -- resize the schematic to that target WIDTH
                # (preserving its own aspect ratio), not a height factor
                # like before, then pad whichever of the two images ends up
                # shorter (their native aspect ratios rarely match) up to
                # the taller one's height.
                SCHEMATIC_WIDTH_FRAC = 0.25 / 0.75  # schematic:row width ratio (1:3)
                target_w = round(row_img.shape[1] * SCHEMATIC_WIDTH_FRAC)
                target_h = round(schematic.shape[0] * target_w / schematic.shape[1])
                schematic = _resize_image(schematic, target_h)
                if schematic.shape[2] != row_img.shape[2]:
                    # Match alpha-channel presence between the two before concatenating.
                    if row_img.shape[2] == 3:
                        row_img = np.concatenate(
                            [row_img, np.ones((*row_img.shape[:2], 1), dtype=row_img.dtype)], axis=2)
                    if schematic.shape[2] == 3:
                        schematic = np.concatenate(
                            [schematic, np.ones((*schematic.shape[:2], 1), dtype=schematic.dtype)], axis=2)
                if row_img.shape[0] < schematic.shape[0]:
                    pad_total = schematic.shape[0] - row_img.shape[0]
                    pad_top = pad_total // 2
                    pad_bottom = pad_total - pad_top
                    row_img = np.concatenate([
                        np.ones((pad_top, row_img.shape[1], row_img.shape[2]), dtype=row_img.dtype),
                        row_img,
                        np.ones((pad_bottom, row_img.shape[1], row_img.shape[2]), dtype=row_img.dtype),
                    ], axis=0)
                elif schematic.shape[0] < row_img.shape[0]:
                    pad_total = row_img.shape[0] - schematic.shape[0]
                    pad_top = pad_total // 2
                    pad_bottom = pad_total - pad_top
                    schematic = np.concatenate([
                        np.ones((pad_top, schematic.shape[1], schematic.shape[2]), dtype=schematic.dtype),
                        schematic,
                        np.ones((pad_bottom, schematic.shape[1], schematic.shape[2]), dtype=schematic.dtype),
                    ], axis=0)
                gap = np.ones((row_img.shape[0], SCHEMATIC_GAP_PX, row_img.shape[2]), dtype=row_img.dtype)
                combined = np.concatenate([schematic, gap, row_img], axis=1)
                mpimg.imsave(output_png, combined)
                log(f"Prepended domain schematic (laser_z={laser_z*1e3:.3f}mm, "
                    f"{time_label} us) to composite")
        else:
            log(f"No schematic script found at {SCHEMATIC_SCRIPT} -- skipping")

        # No separate color->x-position legend anymore (removed, user
        # request, 2026-08-03) -- redundant now that the in-scene
        # "x=..." text labels next to each ghost frame carry the same
        # information directly.
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
    zmin_crop, zmax_crop = Z_VIEW_MIN, Z_VIEW_MAX
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

    # Cell locator on the native mesh, to seed each ray's start-of-line phase
    # via proper interpolation (see start_state() below) rather than
    # snapping to whichever mesh vertex happens to be nearest.
    merged_data = servermanager.Fetch(merged)
    gm_start_arr = vtk_to_numpy(merged_data.GetPointData().GetArray(FIELD_GM))
    ls_start_arr = vtk_to_numpy(merged_data.GetPointData().GetArray(FIELD_LS))
    start_cell_locator = vtk.vtkStaticCellLocator()
    start_cell_locator.SetDataSet(merged_data)
    start_cell_locator.BuildLocator()
    # Fallback for when FindCell misses -- happens more than a rare float
    # edge case (some tens of times per frame, always at a handful of
    # recurring y values -- likely OpenFOAM's general polyhedral cells
    # occasionally tripping up vtkStaticCellLocator's point-in-cell test,
    # or genuine coincidental alignment with an AMR transition boundary
    # repeating at regular z-intervals). Not a regression risk either way:
    # this is exactly the old nearest-vertex behavior, so a ray that hits
    # this path is no worse off than before the fix, just not improved.
    start_point_locator = vtk.vtkPointLocator()
    start_point_locator.SetDataSet(merged_data)
    start_point_locator.BuildLocator()
    log("start-state cell locator built")
    _start_state_fallback_count = [0]

    def start_state(x0, y0, z0):
        """Phase (metal, liquid) at the ray's own start point, via proper
        cell interpolation -- not vtkPointLocator.FindClosestPoint's
        nearest-*vertex* snap, which silently picks the wrong side of the
        interface wherever the local mesh is coarse relative to how close
        the query point sits to the true boundary (visible as a blocky,
        cell-sized misclassification right where the interface passes near
        the fixed x0/x1 slit edge -- user-reported, 2026-08-02). Finding the
        actual enclosing cell and blending its own corner values is the
        standard resolution-correct way to evaluate a field at an arbitrary
        point, and stays entirely local to (x0,y0,z0) -- unlike the
        crossing-parity approach tried and reverted earlier (see task.md),
        this makes no assumption about what phase lies far away, so it
        can't regress the (common, away-from-the-interface) case where the
        solid substrate spans the full domain width with no crossing at
        all.
        """
        cell_id = start_cell_locator.FindCell((x0, y0, z0))
        if cell_id < 0:
            pid = start_point_locator.FindClosestPoint((x0, y0, z0))
            _start_state_fallback_count[0] += 1
            metal = gm_start_arr[pid] >= ISO_THRESHOLD
            liquid = metal and (ls_start_arr[pid] >= LS_TSOLIDUS)
            return metal, liquid
        cell = merged_data.GetCell(cell_id)
        n = cell.GetNumberOfPoints()
        pcoords = [0.0, 0.0, 0.0]
        weights = [0.0] * n
        closest = [0.0, 0.0, 0.0]
        sub_id = vtk.mutable(0)
        dist2 = vtk.mutable(0.0)
        cell.EvaluatePosition((x0, y0, z0), closest, sub_id, pcoords, dist2, weights)
        gm_val = sum(w * gm_start_arr[cell.GetPointId(i)] for i, w in enumerate(weights))
        ls_val = sum(w * ls_start_arr[cell.GetPointId(i)] for i, w in enumerate(weights))
        metal = gm_val >= ISO_THRESHOLD
        liquid = metal and (ls_val >= LS_TSOLIDUS)  # T >= TSolidus
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
    log(f"start_state: FindCell fallback (nearest-vertex) used "
        f"{_start_state_fallback_count[0]} times across {NY * NZ * N_SUB * N_SUB} sub-rays")
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
    # Pin the view to the image's own extent before adding any further
    # overlays below (rays in particular span a much wider x/y/z range in
    # 3D than this crop, and are only reduced to the visible z/y range by
    # matplotlib's axes clipping -- not by pre-filtering the data -- so
    # without this, autoscale could otherwise grow the figure to fit them).
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    # Laser rays (multi-reflection ray-tracing absorption model), if this
    # case has them -- drawn as orange line segments, projected onto this
    # view's (z, y) plane (x dropped, same "collapse the ray-tracing axis"
    # treatment as the attenuation image itself). Each segment's alpha comes
    # from its own endpoints' power (normalized against this frame's own
    # peak), so a ray visibly fades out as it's absorbed/loses energy along
    # its path, rather than being drawn at constant brightness end-to-end.
    # RAY_MAX_OPACITY caps the brightest (freshest, highest-power) segments
    # at 50% rather than fully opaque -- 100% was too bright/dominant over
    # the attenuation image underneath (user feedback, 2026-08-02).
    RAY_MAX_OPACITY = 0.5
    rays = _load_laser_rays(os.path.dirname(foam_file), time_value)
    if rays is not None:
        points, segments, power, ray_idx, rays_vtk_path = rays
        log(f"Loaded laser rays: {rays_vtk_path} ({len(points)} points, "
            f"{len(segments)} segments)")
        if len(segments) and power is not None:
            # Each ray's *recorded* path starts at a fixed launch plane the
            # solver's ray-tracing model uses internally (observed at
            # y~0.2mm for this case) -- not the domain's actual top
            # boundary (ymin, here y=0) and not the true metal surface
            # either. Since that launch height sits inside this view's
            # visible y-window, rays otherwise appear to originate mid-air
            # partway down the frame instead of entering from above
            # (user-reported, 2026-08-02). Add one synthetic segment per
            # ray from (same x/z, ymin) to that first recorded point, at
            # constant (launch) power, so each ray visibly continues up to
            # the domain's own top boundary -- matplotlib's axes clipping
            # (see ax.set_xlim/set_ylim above) takes care of cutting it off
            # at this view's own crop, exactly like every other overlay.
            if ray_idx is not None:
                _, first_idx = np.unique(ray_idx, return_index=True)
                n_orig = len(points)
                launch_points = points[first_idx].copy()
                launch_points[:, 1] = ymin
                launch_power = power[first_idx]
                points = np.concatenate([points, launch_points], axis=0)
                power = np.concatenate([power, launch_power], axis=0)
                launch_segments = np.stack(
                    [np.arange(n_orig, n_orig + len(first_idx)), first_idx], axis=1)
                segments = np.concatenate([segments, launch_segments], axis=0)

            from matplotlib.collections import LineCollection
            p0, p1 = points[segments[:, 0]], points[segments[:, 1]]
            seg_xy = np.stack([
                np.stack([p0[:, 2] * 1e3, p0[:, 1] * 1e3], axis=1),
                np.stack([p1[:, 2] * 1e3, p1[:, 1] * 1e3], axis=1),
            ], axis=1)
            seg_power = (power[segments[:, 0]] + power[segments[:, 1]]) / 2.0
            power_max = power.max()
            alpha = np.clip(seg_power / power_max, 0.0, 1.0) if power_max > 0 else np.zeros_like(seg_power)
            alpha *= RAY_MAX_OPACITY
            colors = np.zeros((len(segments), 4))
            colors[:, :3] = [1.0, 0.55, 0.0]  # orange
            colors[:, 3] = alpha
            ax.add_collection(LineCollection(seg_xy, colors=colors, linewidths=0.5, zorder=1.5))
    else:
        log("No laser-ray VTK series found for this case -- skipping ray overlay")

    ax.plot(zs * 1e3, liquid_bottom * 1e3, linestyle=':', color='skyblue', linewidth=1.2,
            label='melt pool bottom (liquid/solid)', zorder=2)
    ax.axhline(0.4, color='gray', alpha=0.3, linewidth=0.8, zorder=1)  # faint
                                                                        # depth
                                                                        # reference

    # Scale bar overlaid on the image itself (white, so it reads against
    # the solid-black bottom section) instead of numeric z-axis tick
    # labels -- user request, 2026-08-02. Fixed at z=[1.5,2.0]mm (0.5mm),
    # labeled with its own true length; positioned low enough (90% down
    # the y-window) to sit safely below the faint 0.4mm depth-reference
    # line above, comfortably inside the solid region for any frame.
    BAR_FONTSIZE = round(20 * 0.75)  # was a flat 20 -- shrunk to 0.75x
                                       # (user request, 2026-08-02, "the
                                       # size of the bar['s text] to 0.75x
                                       # of the current font size")
    BAR_Z0_MM, BAR_Z1_MM = 1.5, 2.0
    BAR_Y_MM = (Y_DEPTH_MIN + 0.9 * (XRAY_Y_DEPTH_MAX - Y_DEPTH_MIN)) * 1e3
    # Nudged down half a character's height (user request, 2026-08-02),
    # computed from the actual px/mm scale (px_per_mm_y = FIG_HEIGHT_IN*DPI
    # / y_range_mm) rather than a hardcoded mm offset, so it stays correct
    # if the label's fontsize or the figure's DPI/size ever change.
    px_per_mm_y = FIG_HEIGHT_IN * DPI / y_range_mm
    BAR_Y_MM += 0.5 * (BAR_FONTSIZE * DPI / 72) / px_per_mm_y
    CAP_HEIGHT_MM = 0.03 * (XRAY_Y_DEPTH_MAX - Y_DEPTH_MIN) * 1e3
    ax.plot([BAR_Z0_MM, BAR_Z1_MM], [BAR_Y_MM, BAR_Y_MM], color='white',
            linewidth=3, solid_capstyle='butt', zorder=5)
    for zc in (BAR_Z0_MM, BAR_Z1_MM):
        ax.plot([zc, zc], [BAR_Y_MM - CAP_HEIGHT_MM / 2, BAR_Y_MM + CAP_HEIGHT_MM / 2],
                color='white', linewidth=3, zorder=5)
    ax.text(BAR_Z0_MM - 0.04, BAR_Y_MM, f'{BAR_Z1_MM - BAR_Z0_MM:g}mm', color='white',
            fontsize=BAR_FONTSIZE, fontweight='bold', ha='right', va='center', zorder=5)

    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', left=False, labelleft=False)
    # No x-axis ticks/numbers ("z - coord (mm)") -- removed per request,
    # replaced by the scale bar above. bbox_inches='tight' below reclaims
    # the space the axis used to reserve automatically.
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
