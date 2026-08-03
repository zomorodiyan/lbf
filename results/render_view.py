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
# _trim_side_whitespace, _trim_vertical_whitespace, _draw_offset_markers)
# the other three do.
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
GM_RIM_COLORS = [              # one flat color per offset, cyan/green/
                               # yellow (was red/orange/yellow -- user
                               # request, 2026-08-02), indexed in
                               # OFFSETS_BEHIND_LASER order
    [0.0, 1.0, 1.0],             # cyan (was red)
    [0.0, 0.8, 0.0],             # green (was orange; pistachio green
                                  # first tried here didn't read well
                                  # enough against the render -- user
                                  # feedback, 2026-08-02 -- switched to a
                                  # more saturated true green)
    [1.0, 1.0, 0.0],             # yellow (unchanged)
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
                       width_frac=0.35, cb_height=None, bottom_margin=15,
                       right_margin=15, bottom_margin_frac=None,
                       down_shift_chars=0):
    """Render ctf's colorbar *narrower* than the view (width_frac of its
    own already-trimmed width -- deliberately small, not spanning the
    image) and alpha-composite it onto the view image itself, near the
    bottom-right (user request, 2026-08-02, was bottom-center), overlapping
    the actual content there -- not appended
    below (which would grow the canvas) and not a separate legend file
    included via some other composition step. Call *after*
    _trim_side_whitespace()/_trim_vertical_whitespace() (and, for
    render_transverse, after the schematic prepend) so output_png's width
    is already final.

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

    bottom_margin_frac (fraction of vh, overrides bottom_margin when given):
    render_transverse uses this instead of a fixed pixel bottom_margin --
    _render_stacked_video.sh crops the bottom 10% off the transverse row for
    the *stacked composite only* (not this standalone file), so a fixed
    pixel margin would either clear that crop line by a different (and
    possibly too-small) amount depending on this frame's own vh, or get cut
    into outright. A margin expressed as a fraction of vh clears that 10%
    crop by the same proportion on every frame regardless of vh (user
    request, 2026-08-02: "move the colorbar of transverse images enough to
    compensate").

    down_shift_chars: nudge the whole overlaid title+bar block down by this
    many characters (screen-space, same "0.6 * fontsize" character-height
    convention used elsewhere in this repo, e.g. plot_domain_schematic.py's
    axis-label nudges) -- positive moves it down, i.e. reduces its own
    clearance from the bottom edge. User request, 2026-08-02: "move the
    transverse colorbar image 2.5 characters down."
    """
    view_img = mpimg.imread(output_png)
    if view_img.shape[2] == 3:
        alpha = np.ones((*view_img.shape[:2], 1), dtype=view_img.dtype)
        view_img = np.concatenate([view_img, alpha], axis=2)
    vh, vw = view_img.shape[:2]
    if bottom_margin_frac is not None:
        bottom_margin = round(vh * bottom_margin_frac)
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

    # Overlay: standard "over" alpha compositing onto the bottom-right of
    # the view image, in place -- the combined title+bar image's own
    # transparent margins let the underlying view content show through
    # untouched everywhere else.
    cb_img = mpimg.imread(colorbar_png)
    x0 = max(0, vw - right_margin - cb_img.shape[1])
    char_px = 0.6 * label_font_size
    y0 = max(0, vh - bottom_margin - cb_img.shape[0] + round(down_shift_chars * char_px))
    _alpha_paste(view_img, cb_img, x0, y0)
    mpimg.imsave(output_png, view_img)
    log(f"Overlaid colorbar onto: {output_png}")


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
    bounds = _content_col_bounds(img)
    if bounds:
        cmin, cmax = bounds
        mpimg.imsave(output_png, img[:, cmin:cmax + 1])
        log(f"Trimmed side whitespace: kept columns [{cmin},{cmax}] of {img.shape[1]}")
    else:
        log("No content found for side-trim; skipping")


def _content_col_bounds(img, min_col_height_px=5, pad=15):
    """Bounds-only half of _trim_side_whitespace's logic (see its own
    docstring for why a minimum-column-content-height check is needed
    instead of a naive non-white test) -- returns (cmin, cmax) or None,
    without touching any file. Split out so
    _trim_side_whitespace_excluding_markers can compute bounds from one
    rendered image and apply them to a different one."""
    non_white = np.any(img[:, :, :3] < 0.98, axis=2)
    content_cols = np.count_nonzero(non_white, axis=0) >= min_col_height_px
    if not content_cols.any():
        return None
    cmin, cmax = np.where(content_cols)[0][[0, -1]]
    cmin = max(0, cmin - pad)
    cmax = min(img.shape[1] - 1, cmax + pad)
    return cmin, cmax


def _trim_side_whitespace_excluding_markers(view, output_png, marker_disps):
    """Like _trim_side_whitespace, but the crop bounds are computed from a
    version of the scene with the transverse-cut marker lines hidden, so
    those full-crop-height lines don't pull the kept-column range out to
    wherever they land -- e.g. past the real track's own extent early in a
    scan (user-reported, 2026-08-02: "affecting the range that is viewed
    ... not desirable"). Renders twice: once (markers hidden) purely to
    measure bounds, once (markers shown, the real output) to actually
    save -- output_png always ends up as the marker-inclusive render, just
    cropped using the marker-free bounds instead of its own.
    """
    for d in marker_disps:
        d.Visibility = 0
    Render(view)
    SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)
    bounds = _content_col_bounds(mpimg.imread(output_png))

    for d in marker_disps:
        d.Visibility = 1
    Render(view)
    SaveScreenshot(output_png, view, ImageResolution=view.ViewSize)

    if bounds:
        cmin, cmax = bounds
        img = mpimg.imread(output_png)
        mpimg.imsave(output_png, img[:, cmin:cmax + 1])
        log(f"Trimmed side whitespace (markers excluded from bounds): kept columns [{cmin},{cmax}]")
    else:
        log("No content found for side-trim; skipping")


def _trim_vertical_whitespace(output_png):
    """Crop the *saved image* down to its actual content row range.

    FRAME_MARGIN deliberately zooms the camera out a bit so content doesn't
    touch the frame's top/bottom edge (see each view's own camera comment),
    which leaves a blank margin above and below the real crop-window
    content -- trim it off. Unlike _trim_side_whitespace's column check, no
    minimum-run-length trick is needed here: that trick exists because the
    untouched flat plate's edge-on silhouette is a real (if 1px-tall) line
    spanning the *entire* z-window horizontally, regardless of how much
    track exists yet. There's no vertical equivalent -- the blank margin
    here sits entirely *outside* the geometry's own crop box (Y_DEPTH_MIN/
    MAX or X_LATERAL_MIN/MAX), so it's genuinely empty canvas, not a thin
    sliver of real content to filter out.
    """
    img = mpimg.imread(output_png)
    non_white = np.any(img[:, :, :3] < 0.98, axis=2)
    content_rows = np.any(non_white, axis=1)
    if content_rows.any():
        rmin, rmax = np.where(content_rows)[0][[0, -1]]
        pad = 10
        rmin = max(0, rmin - pad)
        rmax = min(img.shape[0] - 1, rmax + pad)
        mpimg.imsave(output_png, img[rmin:rmax + 1, :])
        log(f"Trimmed vertical whitespace: kept rows [{rmin},{rmax}] of {img.shape[0]}")
    else:
        log("No content found for vertical trim; skipping")


def _clip_top_fraction(output_png, frac):
    """Crop off the top `frac` (0-1) of the *already-trimmed* saved image --
    used by render_lateral to cut the top 10% of gas headspace (user
    request, 2026-08-02), independent of _trim_vertical_whitespace's own
    blank-margin trim above."""
    img = mpimg.imread(output_png)
    cut = round(img.shape[0] * frac)
    mpimg.imsave(output_png, img[cut:, :])
    log(f"Clipped top {frac * 100:.0f}%: kept rows [{cut},{img.shape[0] - 1}]")


def _trim_vertical_whitespace_uniform(image_paths):
    """Like _trim_vertical_whitespace, but computes one shared row range
    across *all* given images and applies it to every one -- for
    render_transverse's panels, which must stay the same height to be
    hstacked together afterward (independent per-panel trimming could crop
    each to a different height, since real melt-pool extent genuinely
    differs cut to cut). Takes the union of content rows across panels, so
    nothing real gets clipped in any single panel.
    """
    imgs = [mpimg.imread(p) for p in image_paths]
    rmin_all, rmax_all = None, None
    for img in imgs:
        non_white = np.any(img[:, :, :3] < 0.98, axis=2)
        content_rows = np.any(non_white, axis=1)
        if not content_rows.any():
            continue
        rmin, rmax = np.where(content_rows)[0][[0, -1]]
        rmin_all = rmin if rmin_all is None else min(rmin_all, rmin)
        rmax_all = rmax if rmax_all is None else max(rmax_all, rmax)
    if rmin_all is None:
        log("No content found for vertical trim across panels; skipping")
        return
    pad = 10
    rmin_all = max(0, rmin_all - pad)
    rmax_all = min(imgs[0].shape[0] - 1, rmax_all + pad)
    for p, img in zip(image_paths, imgs):
        mpimg.imsave(p, img[rmin_all:rmax_all + 1, :])
    log(f"Trimmed vertical whitespace (union across panels): kept rows [{rmin_all},{rmax_all}]")


def _trim_side_whitespace_uniform(image_paths, min_col_height_px=5):
    """Like _trim_side_whitespace, but computes one shared column range
    across *all* given images and applies it to every one -- same "stay the
    same size, take the union of content" reasoning as
    _trim_vertical_whitespace_uniform above, just on the other axis. Needed
    because each transverse panel's own camera framing (FRAME_MARGIN) leaves
    blank margin on both sides that a small PANEL_GAP_PX between panels
    can't remove on its own -- panels visibly weren't touching even at a
    tiny gap (user-reported, 2026-08-02). Uses the same minimum-column-
    content-height trick as _trim_side_whitespace (not just "any non-white
    pixel") for the same reason: the untouched flat cap along each panel's
    bottom is a real, non-white row spanning every column already.
    """
    imgs = [mpimg.imread(p) for p in image_paths]
    cmin_all, cmax_all = None, None
    for img in imgs:
        non_white = np.any(img[:, :, :3] < 0.98, axis=2)
        content_cols = np.count_nonzero(non_white, axis=0) >= min_col_height_px
        if not content_cols.any():
            continue
        cmin, cmax = np.where(content_cols)[0][[0, -1]]
        cmin_all = cmin if cmin_all is None else min(cmin_all, cmin)
        cmax_all = cmax if cmax_all is None else max(cmax_all, cmax)
    if cmin_all is None:
        log("No content found for side trim across panels; skipping")
        return
    pad = 6  # narrowed 10 -> 3, then doubled back up to 6 along with
             # PANEL_GAP_PX -- this padding applies on *both* sides of
             # every panel, so it compounds with PANEL_GAP_PX between
             # adjacent panels (user request, 2026-08-02, "the space in
             # between transverse images should twice what it is")
    cmin_all = max(0, cmin_all - pad)
    cmax_all = min(imgs[0].shape[1] - 1, cmax_all + pad)
    for p, img in zip(image_paths, imgs):
        mpimg.imsave(p, img[:, cmin_all:cmax_all + 1])
    log(f"Trimmed side whitespace (union across panels): kept columns [{cmin_all},{cmax_all}]")


def _trim_whitespace_bbox(img, pad=8):
    """Crop `img` (RGB or RGBA array) to the bounding box of its non-white
    content on *all four* sides at once -- unlike _trim_side_whitespace/
    _trim_vertical_whitespace(_uniform) above, which each trim only one
    axis for the ParaView-rendered views (where the other axis is already
    tight or deliberately kept uniform across panels). Used for
    domain_schematic.png, a standalone matplotlib figure with wide margins
    on every side."""
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


def _resize_image(img, new_height):
    """Resize `img` (RGB or RGBA array) to exactly `new_height` px, scaling
    width to preserve its own aspect ratio. Rendered through a throwaway
    matplotlib figure (imshow + savefig-to-buffer) for proper antialiased
    resampling -- same "render through a figure, read back the pixel
    buffer" technique _render_title_image() already uses, since there's no
    PIL/scipy available in this Docker image for a direct resize call."""
    h, w = img.shape[:2]
    new_width = max(1, round(w * new_height / h))
    dpi = 100
    fig = plt.figure(figsize=(new_width / dpi, new_height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, aspect='auto')
    ax.axis('off')
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy().astype(np.float32) / 255.0
    plt.close(fig)
    return buf


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

    Returns the list of marker Show() displays, so callers can temporarily
    hide them (see _trim_side_whitespace_excluding_markers) -- these lines
    span the full crop height by construction, so they trivially satisfy
    _trim_side_whitespace's min-column-content-height check regardless of
    where they land, which otherwise pulls the kept-column range out to
    wherever a marker happens to sit even when it's well past the real
    track's own extent (user-reported, 2026-08-02).
    """
    disps = []
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
        disps.append(disp)
    return disps


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
    marker_disps = _draw_offset_markers(
        view, laser_z,
        lambda z_cut: ([X_LATERAL_MIN, y_marker, z_cut], [X_LATERAL_MAX, y_marker, z_cut]),
    )

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

    _trim_side_whitespace_excluding_markers(view, output_png, marker_disps)
    _trim_vertical_whitespace(output_png)
    _overlay_colorbar(ctf, 'y (mm)', output_png, custom_labels=[-0.1, 0.0, 0.1])

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
    log(f"Mushy-zone (liquid/solid) surface: {ls_poly.GetNumberOfCells()} cells")

    # Exact plane intersection at x=0 -- not a rendered clip, so there's no
    # depth-buffer ambiguity to arbitrate for this 1D curve (same technique
    # as render_transverse's slice_at_cut()).
    ls_slice = Slice(Input=ls_contour)
    ls_slice.SliceType = 'Plane'
    ls_slice.SliceType.Origin = [0.0, (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0, (z_window_min + z_window_max) / 2.0]
    ls_slice.SliceType.Normal = [1.0, 0.0, 0.0]
    ls_slice.UpdatePipeline(time=time_value)

    # Bug fix: ls_slice comes from ls_contour, which is built off the *whole*
    # domain (unclipped, unlike gm_contour -> feature below) -- the x=0
    # slice through it can carry mushy-zone content at any z/y across the
    # entire mesh, not just inside [z_window_min,z_window_max] x
    # [Y_DEPTH_MIN,Y_DEPTH_MAX]. Camera framing alone doesn't stop this from
    # showing up: FRAME_MARGIN's zoom-out margin means the actual rendered
    # z-range is wider than [z_window_min,z_window_max] (by the same
    # factor), so real mushy-outline content just past the intended window
    # was visible in the saved image, undermining the whole point of fixing
    # Z_VIEW_MIN/MAX. Explicit box-clip, matching feature's own bounds.
    ls_slice_clip = Clip(Input=ls_slice)
    ls_slice_clip.ClipType = 'Box'
    ls_slice_clip.ClipType.Position = [-1e-6, Y_DEPTH_MIN, z_window_min]
    ls_slice_clip.ClipType.Length = [2e-6, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
    ls_slice_clip.Invert = 1
    ls_slice_clip.UpdatePipeline(time=time_value)
    ls_slice_poly = servermanager.Fetch(ls_slice_clip)
    log(f"Mushy-zone outline at x=0 (clipped to view window): {ls_slice_poly.GetNumberOfCells()} cells "
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

    marker_disps = _draw_offset_markers(
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

    _trim_side_whitespace_excluding_markers(view, output_png, marker_disps)
    _trim_vertical_whitespace(output_png)
    _clip_top_fraction(output_png, 0.10)
    _overlay_colorbar(ctf, 'x (mm)', output_png, custom_labels=[-0.1, 0.0, 0.1])

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
    CAP_SOLID_COLOR = [0.9, 0.9, 0.9]      # flat light gray, filled solid
                                            # portion of the cross-section
                                            # cap -- lightened from 0.8
                                            # (user request, 2026-08-02)
    LS_SURFACE_COLOR = [0.75, 0.75, 0.75]  # flat lighter gray, mushy-zone
                                            # surface itself (not colored)
    PANEL_GAP_PX = 6                        # white gap between panels --
                                            # narrowed 20 -> 6 -> 3, then
                                            # doubled back up to 6 (user
                                            # request, 2026-08-02, "the
                                            # space in between transverse
                                            # images should twice what it
                                            # is") -- the per-panel side
                                            # trim pad below is doubled to
                                            # match, so the total visible
                                            # gap (2*pad + this) scales by
                                            # the same factor
    SCHEMATIC_GAP_PX = 0                   # white gap between the schematic
                                            # and the first panel -- narrowed
                                            # 20 -> 8 -> 4 -> 0 (user request,
                                            # 2026-08-02, "no need to have
                                            # any gap in between transverse
                                            # images and the schematic")
    SCHEMATIC_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     'plot_domain_schematic.py')  # regenerated
                                     # per-frame below (laser position/time
                                     # label vary), not a static pre-baked
                                     # image
    # Local override, transverse-only (shadows the shared module-level
    # Y_DEPTH_MIN/MAX used by top/lateral, same "give this view its own
    # value" pattern as render_xray's XRAY_Y_DEPTH_MAX) -- user request,
    # 2026-08-02: each panel should show 100-450um below the surface (was
    # 50-200um, then 100-500um).
    Y_DEPTH_MIN = 0.1e-3
    Y_DEPTH_MAX = 0.45e-3

    # Color-scale bounds for z_minus_laser = (z - laser_z)*1e3 (mm) -- always
    # <= 0 within the cropped region (z is clamped to <= z_cut < laser_z).
    # Z_REL_MIN is the most-negative end (furthest offset + 1mm further out,
    # for headroom), Z_REL_MAX is the least-negative end (the first-listed/
    # largest offset, i.e. closest to 0) -- the whole panel set shares one
    # consistent, comparable scale.
    Z_REL_MIN = -(max(OFFSETS_BEHIND_LASER) * 1e3 + 1.0)
    Z_REL_MAX = -(OFFSETS_BEHIND_LASER[0] * 1e3)

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
            gm_rim_color = GM_RIM_COLORS[panel_idx % len(GM_RIM_COLORS)]
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

            # z_minus_laser (mm) for coloring the gas/metal surface --
            # each point's own z minus laser_z (<= 0 throughout this crop).
            gm_zcolor = Calculator(Input=feature)
            gm_zcolor.AttributeType = 'Point Data'
            gm_zcolor.ResultArrayName = 'z_minus_laser'
            gm_zcolor.Function = f'(coordsZ-{laser_z})*1e3'
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
            ColorBy(disp, ('POINTS', 'z_minus_laser'))

            ctf = GetColorTransferFunction('z_minus_laser')
            # Smooth white->blue gradient across the full Z_REL_MIN-MAX
            # range (ascending X order: most-negative/furthest-from-laser
            # first). White is the Z_REL_MIN endpoint, not a midpoint on the
            # way to red -- ParaView's default diverging interpolation would
            # otherwise wash out through white *then* continue past it.
            ctf.RGBPoints = [
                Z_REL_MIN, 1.0, 1.0, 1.0,
                Z_REL_MAX, 0.0, 0.0, 1.0,
            ]

            # Cross-section cap's solid portion, drawn last (frontmost --
            # it sits exactly at the clip window's own front face, so this
            # only matters for exact ties, not real occlusion). Back to a
            # flat gray fill (was briefly swapped to the panel's own
            # gm_rim_color, then reverted -- user request, 2026-08-02,
            # "color the solid surface as gray like before"); the
            # crop-boundary side outlines below stay white.
            cap_solid_disp = Show(cap_solid, view)
            cap_solid_disp.Representation = 'Surface'
            cap_solid_disp.ColorArrayName = ['POINTS', '']
            cap_solid_disp.AmbientColor = CAP_SOLID_COLOR
            cap_solid_disp.DiffuseColor = CAP_SOLID_COLOR

            # Straight sides of the solid cap (the box-clip crop boundary --
            # left/right/bottom of the crop window, as opposed to the
            # actual curvy gas/metal interface at the top, already outlined
            # by gm_rim below) colored plain white -- now that the cap fill
            # itself carries the panel's cross-section color, the crop-
            # boundary outline is white to stay visually distinct from it
            # (swapped from the panel's own rim color -- user request,
            # 2026-08-02). 3 explicit Line sources (same pattern as the
            # other markers in this function) rather than enabling
            # EdgeVisibility on cap_solid_disp itself, which would draw
            # every individual mesh cell's edges (a busy internal
            # wireframe), not just this outer silhouette.
            for p1, p2 in [
                ((X_LATERAL_MIN, Y_DEPTH_MIN, z_cut), (X_LATERAL_MIN, Y_DEPTH_MAX, z_cut)),  # left
                ((X_LATERAL_MAX, Y_DEPTH_MIN, z_cut), (X_LATERAL_MAX, Y_DEPTH_MAX, z_cut)),  # right
                ((X_LATERAL_MIN, Y_DEPTH_MAX, z_cut), (X_LATERAL_MAX, Y_DEPTH_MAX, z_cut)),  # bottom
            ]:
                cap_side = Line(Point1=p1, Point2=p2)
                cap_side_disp = Show(cap_side, view)
                cap_side_disp.Representation = 'Wireframe'
                cap_side_disp.ColorArrayName = ['POINTS', '']
                cap_side_disp.AmbientColor = [1.0, 1.0, 1.0]
                cap_side_disp.DiffuseColor = [1.0, 1.0, 1.0]
                cap_side_disp.LineWidth = GM_RIM_LINE_WIDTH

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

        # Trim blank top/bottom margin (from FRAME_MARGIN's camera zoom-out,
        # same mechanism as the other two views) -- union across all 3
        # panels so they stay the same height for the hstack below.
        _trim_vertical_whitespace_uniform(panel_pngs)
        # Same idea, sideways: each panel's own camera framing leaves blank
        # margin left/right too, which PANEL_GAP_PX alone can't remove --
        # panels weren't actually touching even at a small gap (user-
        # reported, 2026-08-02). Union across panels, same as above, so
        # they stay the same (now narrower) width for the hstack below.
        _trim_side_whitespace_uniform(panel_pngs)

        # Stack the panels side-by-side -- now genuinely tight against
        # PANEL_GAP_PX, not just nominally the same size.
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

        # Prepend the domain schematic to the left of the 3 cut panels
        # *before* the colorbar overlay below -- unlike the previous
        # ordering (colorbar first, schematic prepended after), the
        # colorbar's own size now needs to be computed against the row's
        # true final width (see the comment on _overlay_colorbar's call
        # below), which only exists once the schematic is already in
        # place. Regenerated fresh for *this* frame (not a static pre-baked
        # image) -- its
        # single laser marker/time label must match the actual timestep
        # being rendered here, so plot_domain_schematic.py is invoked as a
        # subprocess with this frame's own laser_z and a time label, output
        # to tmp_dir (cleaned up automatically below). Cropped to its own
        # content bounding box (its source figure has wide margins on every
        # side) and resized to the row's own height, preserving aspect ratio.
        if os.path.exists(SCHEMATIC_SCRIPT):
            schematic_png = os.path.join(tmp_dir, 'schematic.png')
            time_label = f"{time_value * 1e6:.0f} us"
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
                row_img = mpimg.imread(output_png)
                schematic = mpimg.imread(schematic_png)
                schematic = _trim_whitespace_bbox(schematic)
                # Resized noticeably taller than the panel row itself (was
                # matched 1:1) -- at 1:1 the schematic read too small/hard
                # to make out (user report, 2026-08-02), so it's now the
                # larger element and the panel row gets padded (white,
                # top/bottom, centered) up to match its height instead of
                # the other way around.
                SCHEMATIC_HEIGHT_FACTOR = 1.35
                schematic = _resize_image(schematic, round(row_img.shape[0] * SCHEMATIC_HEIGHT_FACTOR))
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
                gap = np.ones((row_img.shape[0], SCHEMATIC_GAP_PX, row_img.shape[2]), dtype=row_img.dtype)
                combined = np.concatenate([schematic, gap, row_img], axis=1)
                mpimg.imsave(output_png, combined)
                log(f"Prepended domain schematic (laser_z={laser_z*1e3:.3f}mm, "
                    f"{time_label}) to composite")
        else:
            log(f"No schematic script found at {SCHEMATIC_SCRIPT} -- skipping")

        # One shared colorbar for z_minus_laser -- same color scale
        # (Z_REL_MIN..Z_REL_MAX) across all 3 panels, so one overlay on the
        # final row is enough; `ctf` here is whichever panel's loop
        # iteration ran last, but GetColorTransferFunction('z_minus_laser')
        # returns the same shared transfer-function proxy every time
        # regardless. 3 explicit labels (min/mid/max), same "readable at a
        # glance" convention as top/lateral's fixed +/-100um -- ParaView's
        # auto-picked labels here were a cluttered 8+ ticks (Z_REL_MIN/MAX
        # aren't round numbers). Overlaid *after* the schematic prepend
        # (was before) -- _overlay_colorbar sizes itself proportionally to
        # output_png's own width so it comes out a consistent size once
        # _render_stacked_video.sh later scales every view to one shared
        # width (user report, 2026-08-02); that only works if the width it
        # measures here is the row's true final width, schematic included,
        # matching what top/lateral/xray already are when their own
        # colorbars get overlaid. bottom_margin_frac (fraction of vh, not a
        # fixed pixel bottom_margin) -- transverse's panels are trimmed
        # tight with no blank margin at the bottom (unlike top/lateral), so
        # a small fixed pixel margin left the bar sitting flush against the
        # crop window's own flat bottom edge, reading as clipped (user
        # report, 2026-08-02). 0.28 (was 0.18) clears
        # _render_stacked_video.sh's own bottom-20%-of-the-row crop
        # (10% -> 20%, applied only to the stacked composite, not this
        # file) by the same ~0.08 comfortable margin as before, regardless
        # of this frame's own vh (user report, 2026-08-02: "move the
        # colorbar ... enough to compensate" / "more cut-off from bottom
        # of the transverse cross section views" -- keep this in sync with
        # that script's own crop fraction if either changes).
        _overlay_colorbar(
            ctf, 'z - z_l (mm)', output_png,
            custom_labels=[round(Z_REL_MIN, 1), round((Z_REL_MIN + Z_REL_MAX) / 2, 1), round(Z_REL_MAX, 1)],
            bottom_margin_frac=0.28, down_shift_chars=2.5,
        )
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
