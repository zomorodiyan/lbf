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
import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import numpy as np

# Liberation Sans -- a free (SIL Open Font License), metrically-compatible
# substitute for Arial (user request, 2026-08-22: "use a font very similar
# to arial if using arial is illegal otherwise download it" -- true Arial
# is a proprietary Monotype font with no legal free-redistribution license,
# so this is the standard legal substitute instead, not Arial itself).
# This Docker image's bundled matplotlib only ships DejaVu/STIX/Computer
# Modern, so the actual font files are downloaded once (see
# results/fonts/README, if present) and registered here at import time
# rather than relying on any system font install.
_font_paths = [os.path.join('/workspace/results/fonts', _f) for _f in (
    'LiberationSans-Regular.ttf', 'LiberationSans-Bold.ttf',
    'LiberationSans-Italic.ttf', 'LiberationSans-BoldItalic.ttf')]
# addfont() doesn't exist yet in this image's matplotlib (3.1.1, predates
# 3.2's addfont) -- createFontList()+ttflist.extend() is that version's
# own equivalent for registering fonts outside the standard search paths.
fm.fontManager.ttflist.extend(fm.createFontList(_font_paths))
plt.rcParams['font.family'] = 'Liberation Sans'

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
CROSS_SECTION_Y_WINDOW = (0.225e-3, 0.375e-3)  # m -- bottom (deeper) edge extended
                               # 20% taller than the original 0.35e-3, top
                               # (shallow) edge unchanged (user request,
                               # 2026-08-11). Matches render_transverse's own
                               # Y_MIN/Y_MAX -- keep in sync.
CROSS_SECTION_MARKER_COLORS = [                # dark blue/green/red, matches
    [0.051, 0.278, 0.631],                       # render_transverse's own
    [0.145, 0.392, 0.157],                       # DARK_COLORS hexes
    [0.718, 0.110, 0.110],                       # (#0d47a1/#256428/#b71c1c)
]
GM_RIM_LINE_WIDTH = 3.0       # wireframe line width (px) for those markers
CBAR_MARGIN_MM = 0.025        # render_cutaway (mm): how far past its own
                               # gradient transition bounds its colorbar's
                               # *visual* range extends -- kept separate
                               # from the much wider domain-derived range
                               # used for actual data coloring, so the bar
                               # itself isn't mostly flat blue/red with the
                               # real gradient squeezed into a sliver (user
                               # request, 2026-08-17). render_top uses its
                               # own local CBAR_MARGIN_Y_MM (+-50um) instead
                               # -- diverged from this one, same date.
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
    if not series.get('files'):
        # Stale series file (empty "files" list) -- happens after any
        # pause/resume until tutorials/laserbeamFoam/fix_vtk_series.py is
        # rerun (see CLAUDE.md). Same "no ray data available" outcome as the
        # file not existing at all -- degrade gracefully rather than crash
        # the whole render on min() over an empty sequence.
        return None
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
                       down_shift_chars=0, side='right', vert='top',
                       labels_above=False, skip_overlay=False, thickness_scale=1.0):
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

    labels_above: put the tick value labels above the bar strip instead of
    below (the default, used by every caller before render_cutaway) --
    user request, 2026-08-17, for render_cutaway specifically. Only
    affects the ScalarBar's own TextPosition; the title stays to the left
    of the bar either way (see this function's own note on why above).

    skip_overlay: save the standalone title+bar colorbar_png (as always),
    but don't alpha-paste it onto output_png itself -- for render_cutaway
    (user request, 2026-08-17): its colorbar is meant to visually overlap
    into the *adjacent* panel above once stacked, which is impossible to
    do from here (this function only ever sees output_png's own canvas,
    and the panels don't get combined until the separate ffmpeg stacking
    step). The stacking script instead overlays the saved standalone file
    itself, positioned to straddle the seam between the two already-
    stacked panels.
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
    # thickness_scale: render_cutaway-only knob -- user request, 2026-08-22:
    # "twice the current thickness" for its topright colorbar, independent
    # of the collage script's own overall-image-size scaling (150%->160%),
    # which just rescales the whole already-rendered PNG uniformly rather
    # than the bar's own stroke thickness. Only touches ScalarBarThickness
    # (leaving it at ParaView's own default) when explicitly requested --
    # multiplying by font_scale unconditionally would also silently rescale
    # every *other* existing caller's bar thickness in proportion to their
    # own view width (font_scale is ~vw/REF_VW, essentially never exactly
    # 1.0), which nothing asked for.
    if thickness_scale != 1.0:
        colorbar.ScalarBarThickness = round(colorbar.ScalarBarThickness * thickness_scale * font_scale)
    label_font_size = colorbar.LabelFontSize  # captured for down_shift_chars
                                               # below -- cb_view/colorbar
                                               # get Delete()d before then
    # Tick value labels below the bar strip by default, not above it (the
    # ParaView default) -- user request, 2026-08-02. labels_above=True
    # (render_cutaway only, 2026-08-17) flips this back to above.
    colorbar.TextPosition = ('Ticks right/top, annotations left/bottom' if labels_above
                              else 'Ticks left/bottom, annotations right/top')
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

    if skip_overlay:
        log(f"skip_overlay=True: left {output_png} untouched, standalone "
            f"colorbar saved to {colorbar_png} for the caller to composite "
            f"elsewhere")
        return

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


def _add_cutaway_grid(output_png, z_left_mm, z_right_mm, y_top_offset_um, y_bottom_offset_um,
                       z_step_mm=0.1, y_step_um=50, time_us=None):
    """Annotation-planning aid for render_cutaway (user request,
    2026-08-18): overlay a physical-coordinate grid -- z in mm, y in um
    offset from the nominal surface (same offset convention as render_top's
    own y-coloring, positive = below surface) -- on a *copy* of output_png,
    saved alongside it as ..._grid.png, so annotation locations/timesteps
    can be discussed by coordinate instead of eyeballed pixels. Leaves
    output_png itself untouched.

    z_left_mm/z_right_mm/y_top_offset_um/y_bottom_offset_um must be the
    *actual* extent this frame was rendered at (the caller derives them
    from the same CameraParallelScale/crop values used for the real
    render) -- passed in rather than re-derived here so there's exactly
    one place (the caller) that has to get that math right, and this
    function just trusts it. Uses imshow(..., extent=...) so matplotlib's
    own coordinate transform places the gridlines/ticks, instead of hand
    computing pixel positions a second time here.
    """
    img = mpimg.imread(output_png)
    h, w = img.shape[:2]
    dpi = 100
    fig, ax = plt.subplots(figsize=(w / dpi * 1.12, h / dpi * 1.18), dpi=dpi)
    ax.imshow(img, extent=(z_left_mm, z_right_mm, y_bottom_offset_um, y_top_offset_um), aspect='auto')
    z_ticks = np.arange(np.ceil(z_left_mm / z_step_mm) * z_step_mm, z_right_mm, z_step_mm)
    y_ticks = np.arange(np.ceil(y_top_offset_um / y_step_um) * y_step_um, y_bottom_offset_um, y_step_um)
    # Each gridline gets its own color, sampled from a position along a
    # colormap -- not a uniform color for every line -- so a line can be
    # picked out and named by eye ("the orange vertical line") instead of
    # having to read its printed number (user request, 2026-08-18). z
    # (vertical lines) and y (horizontal lines) use different colormap
    # families (rainbow vs cool) so which *axis* a color belongs to is
    # also unambiguous at a glance, not just which position along it.
    z_cmap, y_cmap = plt.get_cmap('rainbow'), plt.get_cmap('cool')
    z_span = (z_right_mm - z_left_mm) or 1.0
    y_span = (y_bottom_offset_um - y_top_offset_um) or 1.0
    z_colors = [z_cmap((v - z_left_mm) / z_span) for v in z_ticks]
    y_colors = [y_cmap((v - y_top_offset_um) / y_span) for v in y_ticks]
    for zt, c in zip(z_ticks, z_colors):
        ax.axvline(zt, color=c, linewidth=1.1, alpha=0.9)
    for yt, c in zip(y_ticks, y_colors):
        ax.axhline(yt, color=c, linewidth=1.1, alpha=0.9)
    ax.set_xticks(z_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels([f"{v:.1f}" for v in z_ticks], fontsize=7)
    ax.set_yticklabels([f"{v:+.0f}" for v in y_ticks], fontsize=7)
    for label, c in zip(ax.get_xticklabels(), z_colors):
        label.set_color(c)
        label.set_fontweight('bold')
    for label, c in zip(ax.get_yticklabels(), y_colors):
        label.set_color(c)
        label.set_fontweight('bold')
    ax.set_xlabel("z (mm)", fontsize=8)
    ax.set_ylabel("y offset from surface (um)", fontsize=8)
    if time_us is not None:
        # Top-left, in axes-fraction coords (user request, 2026-08-18:
        # "integer us" time-progression label).
        ax.text(0.015, 0.97, f"{round(time_us)}us", transform=ax.transAxes,
                fontsize=11, fontweight='bold', color='black', ha='left', va='top')
    grid_png = re.sub(r'\.png$', '_grid.png', output_png)
    fig.tight_layout()
    fig.savefig(grid_png)
    plt.close(fig)
    log(f"Saved grid overlay: {grid_png}")


def _overlay_rays_on_cutaway(output_png, case_dir, time_value, ymin_domain,
                              z_left_mm, z_right_mm, y_top_offset_um, y_bottom_offset_um,
                              supersample=1, ray_max_opacity=0.25, ray_color=(0.55, 0.0, 0.0)):
    """Composite laser ray-tracing segments onto output_png in place, dark
    red (user request, 2026-08-22 -- was dark orange first, but that was
    too close to the surface's own orange x-coloring; was green before
    that), each segment's opacity scaled 0-ray_max_opacity (25%, was 50%,
    same request) by its own power (user request, 2026-08-22, for
    collage_v3) -- reuses the exact projected-segments/power-alpha technique
    render_xray's own (orange) ray overlay already uses (see its comments),
    just composited directly onto this ParaView cutaway render in its
    y-offset-from-surface/z-mm convention instead of a separate matplotlib
    figure in raw-y/mm.

    Unlike _add_cutaway_grid (which deliberately saves a separate, differently
    -sized diagnostic file), this OVERWRITES output_png at its own original
    pixel dimensions (fig.add_axes([0,0,1,1]), no tight_layout/margins) --
    output_png is the actual collage source panel, and downstream collage
    cropping (_build_cutaway_collage_v1.py's build_panel) computes crop
    pixel bounds analytically from the image's own shape, so its dimensions
    must not change here.
    """
    rays = _load_laser_rays(case_dir, time_value)
    if rays is None:
        log("No laser-ray VTK series found for this case -- skipping ray overlay")
        return
    points, segments, power, ray_idx, rays_vtk_path = rays
    log(f"Loaded laser rays for cutaway overlay: {rays_vtk_path} "
        f"({len(points)} points, {len(segments)} segments)")
    if not len(segments) or power is None:
        log("Ray VTK has no segments/power data -- skipping ray overlay")
        return

    SURFACE_Y = 0.2e-3
    # Extend each ray from its first *recorded* point up to the domain's
    # real top edge (ymin_domain) at constant (launch) power, same "ray
    # otherwise appears to start mid-air" fix render_xray already applies
    # (see its own comment) -- this view's own set_xlim/set_ylim below
    # clips it to the visible window exactly the same way.
    if ray_idx is not None:
        _, first_idx = np.unique(ray_idx, return_index=True)
        n_orig = len(points)
        launch_points = points[first_idx].copy()
        launch_points[:, 1] = ymin_domain
        launch_power = power[first_idx]
        points = np.concatenate([points, launch_points], axis=0)
        power = np.concatenate([power, launch_power], axis=0)
        launch_segments = np.stack(
            [np.arange(n_orig, n_orig + len(first_idx)), first_idx], axis=1)
        segments = np.concatenate([segments, launch_segments], axis=0)

    from matplotlib.collections import LineCollection
    p0, p1 = points[segments[:, 0]], points[segments[:, 1]]
    seg_xy = np.stack([
        np.stack([p0[:, 2] * 1e3, (p0[:, 1] - SURFACE_Y) * 1e6], axis=1),
        np.stack([p1[:, 2] * 1e3, (p1[:, 1] - SURFACE_Y) * 1e6], axis=1),
    ], axis=1)
    seg_power = (power[segments[:, 0]] + power[segments[:, 1]]) / 2.0
    power_max = power.max()
    alpha = np.clip(seg_power / power_max, 0.0, 1.0) if power_max > 0 else np.zeros_like(seg_power)
    alpha *= ray_max_opacity
    colors = np.zeros((len(segments), 4))
    colors[:, :3] = ray_color
    colors[:, 3] = alpha

    img = mpimg.imread(output_png)
    h, w = img.shape[:2]
    dpi = 100
    fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(img, extent=(z_left_mm, z_right_mm, y_bottom_offset_um, y_top_offset_um), aspect='auto')
    ax.set_xlim(z_left_mm, z_right_mm)
    ax.set_ylim(y_bottom_offset_um, y_top_offset_um)
    # Linewidth scaled with supersample, same reasoning as render_cutaway's
    # own outline_tube_radius/world_per_px -- this is drawn in matplotlib
    # points, which don't otherwise track a supersampled image's extra
    # pixel density (fig's dpi is fixed at 100 regardless of supersample).
    ax.add_collection(LineCollection(seg_xy, colors=colors, linewidths=0.6 * supersample))
    ax.axis('off')
    fig.savefig(output_png, dpi=dpi)
    plt.close(fig)
    log(f"Overlaid {len(segments)} ray segments (green, up to "
        f"{ray_max_opacity * 100:.0f}% opacity by power)")


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
    """Dashed line marker for CROSS_SECTION_X[0] (the -25um cut the
    cutaway view uses), spanning [z_window_min, z_window_max] at that
    fixed x, placed at y=y_marker -- nearer the camera than any real
    geometry (see render_top's own header), so never occluded.

    Black, dashed (user request, 2026-08-17, replacing the earlier solid
    dark-blue marker -- "change the cross-section blue line to a better
    representative... like a very professional article figure... this
    line is defining the second view"): a dashed reference/cutting-plane
    line is the standard convention in technical and scientific figures
    for "this marks a section plane," distinct from a solid line (which
    reads as a real geometric edge) -- and black stays legible against
    the top view's own blue/purple/orange/red data coloring, where a blue
    line would visually compete with it.

    Built as many short, separate Line segments with gaps, not one
    continuous line with a dash *style* applied -- ParaView's modern
    (OpenGL2-era) rendering backend has no supported way to dash a
    Wireframe line through the Show()/Representation API (line stippling
    was a fixed-function-pipeline feature, dropped in the OpenGL2
    rewrite; matches this file's own experience elsewhere, e.g.
    render_lateral's outline lines are solid-only for the same reason).

    Was 3 markers (one per CROSS_SECTION_X, spanning the narrower
    CROSS_SECTION_Z_WINDOW) -- reduced to just this one, and extended to
    the caller's full z-range instead (user request, 2026-08-16: drop the
    x=0/+25um markers and the old narrow z-window, both leftovers from
    render_transverse's old 3-cut design, which isn't what's being
    cross-referenced here anymore -- extend the remaining line the full
    z-range so it reads as "the plane the cutaway view cuts along").
    render_top now calls this with Z_VIEW_MIN/MAX, not CROSS_SECTION_Z_WINDOW.

    Flat-colored, not scalar-colored: ColorArrayName=['POINTS',''], not
    ColorBy(rep, None) -- the latter crashes when there's no array to
    default to (relevant here since these Line sources carry no data).
    """
    x_cut = CROSS_SECTION_X[0]
    DASH_COLOR = [0.0, 0.0, 0.0]
    DASH_LEN_M = 0.06e-3
    GAP_LEN_M = 0.04e-3
    n_dashes = 0
    z = z_window_min
    while z < z_window_max:
        z_end = min(z + DASH_LEN_M, z_window_max)
        marker = Line(Point1=[x_cut, y_marker, z], Point2=[x_cut, y_marker, z_end])
        disp = Show(marker, view)
        disp.Representation = 'Wireframe'
        disp.ColorArrayName = ['POINTS', '']
        disp.AmbientColor = DASH_COLOR
        disp.DiffuseColor = DASH_COLOR
        disp.LineWidth = GM_RIM_LINE_WIDTH * 2  # 2x the shared constant, just
                                                 # for this marker (user
                                                 # request, 2026-08-17) -- not
                                                 # changed at its source, since
                                                 # GM_RIM_LINE_WIDTH is also
                                                 # used by
                                                 # _draw_cross_section_frame_lateral,
                                                 # which wasn't asked to change
        z += DASH_LEN_M + GAP_LEN_M
        n_dashes += 1
    log(f"Cross-section marker: x={x_cut*1e6:+.1f}um -> {n_dashes} dashes at "
        f"y={y_marker*1e3:.3f}mm, z=[{z_window_min*1e3:.3f},{z_window_max*1e3:.3f}]mm")


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
    # Explicit asymmetric transition bounds -- -100um to +200um (user
    # request, 2026-08-17, replacing the earlier symmetric +-0.15mm
    # attempt) -- close to the actual logged content bounds across several
    # tr69 timesteps (roughly -0.08/-0.09mm spatter to +0.16/+0.17mm
    # depression, see git history for that earlier attempt's own notes),
    # rounded to clean numbers with a bit of margin on the depression side.
    # y=0 (the true surface) is no longer exactly the transition's
    # midpoint, so it won't render as pure white -- it lands about 1/3 of
    # the way from blue to red instead. ymin_off/ymax_off (the true domain
    # bounds) stay as the fully-saturated flat endpoints beyond this band.
    TRANSITION_BLUE_MM = -0.100
    TRANSITION_RED_MM = 0.100
    ctf.RGBPoints = [
        ymin_off,            0.0, 0.0, 1.0,
        TRANSITION_BLUE_MM,  0.0, 0.0, 1.0,
        TRANSITION_RED_MM,   1.0, 0.0, 0.0,
        ymax_off,            1.0, 0.0, 0.0,
    ]

    y_marker = ymin - 0.02 * (ymax - ymin)
    _draw_cross_section_markers_top(view, y_marker, Z_VIEW_MIN, Z_VIEW_MAX)

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
    #
    # The colorbar's own visual range tracks the transition bounds, not
    # the full domain-derived ymin_off/ymax_off `ctf` itself uses for
    # actual data coloring (user request, 2026-08-17) -- otherwise the
    # bar is mostly flat blue/red with the real gradient squeezed into a
    # sliver. A separate CTF (never bound to any real Show()/ColorBy --
    # GetScalarBar just reads its RGBPoints directly) keeps this purely
    # cosmetic without touching how ycolor is actually colored.
    #
    # +-50um margin here, not the shared module-level CBAR_MARGIN_MM
    # (+-25um, still used by render_cutaway's own colorbar) -- diverged
    # per-view, user request, 2026-08-17.
    #
    # Labeled in um, not mm (user request, 2026-08-17) -- cbar_ctf is
    # purely cosmetic (never bound to any real Show()/ColorBy, see above),
    # so scaling its own RGBPoints/labels by 1000x here has no effect on
    # how ycolor itself is actually colored (that still happens through
    # the separate, untouched `ctf`, still in mm).
    UM_PER_MM = 1000.0
    CBAR_MARGIN_Y_MM = 0.050
    cbar_ctf = GetColorTransferFunction(FIELD_COLOR + '_cbar_display')
    cbar_ctf.RGBPoints = [
        (TRANSITION_BLUE_MM - CBAR_MARGIN_Y_MM) * UM_PER_MM, 0.0, 0.0, 1.0,
        TRANSITION_BLUE_MM * UM_PER_MM,                      0.0, 0.0, 1.0,
        TRANSITION_RED_MM * UM_PER_MM,                       1.0, 0.0, 0.0,
        (TRANSITION_RED_MM + CBAR_MARGIN_Y_MM) * UM_PER_MM,  1.0, 0.0, 0.0,
    ]
    _overlay_colorbar(cbar_ctf, 'y (μm)', output_png,
                       custom_labels=[TRANSITION_BLUE_MM * UM_PER_MM, 0.0, TRANSITION_RED_MM * UM_PER_MM])

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
    # Solidus/liquidus contour-line color, fixed across all 3 panels
    # (user request, 2026-08-11: "#222 lines at the contourlines of
    # T_solidus and T_liquidus") -- deliberately *not* the per-cut
    # dark_color used for the gas/metal rim below: that one is hued per
    # cross-section so the 3 panels stay identifiable, but solidus/liquidus
    # are the same physical boundary drawn identically on all 3, so a
    # single fixed dark gray reads as "one consistent overlay" rather than
    # blending into whichever hue that panel's fill happens to use.
    # Solidus solid, liquidus dashed (same dash pattern as the ghost-frame
    # leader line below) so the two stay distinguishable at a glance
    # despite sharing a color.
    MELT_BOUNDARY_COLOR = '#222222'
    MELT_BOUNDARY_LINE_WIDTH = 1.2
    MELT_BOUNDARY_LIQUIDUS_DASH = (0, (3, 2))
    # z window (m), fixed in absolute z (not laser-relative) so cross-
    # sections stay comparable across an entire batch of timesteps.
    Z_WINDOW_MIN, Z_WINDOW_MAX = 1.5e-3, 2.0e-3
    # y window (m), shared across all cross-sections. Bottom (deeper) edge
    # extended 20% taller than the original 0.35e-3, top (shallow) edge
    # left unchanged (user request, 2026-08-11) -- matches the shared
    # CROSS_SECTION_Y_WINDOW module constant above, keep in sync.
    Y_MIN, Y_MAX = 0.225e-3, 0.375e-3

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

    # Solidus/liquidus contours (drawn per-cut below as #222222 overlay
    # lines, user request, 2026-08-11), restricted to inside the metal
    # first -- same Clip-then-Contour pattern render_lateral uses (see its
    # own comment on why Clip, not Threshold).
    metal_only = Clip(Input=merged)
    metal_only.ClipType = None
    metal_only.Scalars = ['POINTS', FIELD_GM]
    metal_only.Value = ISO_THRESHOLD
    metal_only.Invert = 0
    metal_only.UpdatePipeline(time=time_value)

    ls_solidus_contour = Contour(Input=metal_only)
    ls_solidus_contour.ContourBy = ['POINTS', FIELD_LS]
    ls_solidus_contour.Isosurfaces = [LS_TSOLIDUS]
    ls_solidus_contour.UpdatePipeline(time=time_value)
    ls_solidus_poly = servermanager.Fetch(ls_solidus_contour)
    log(f"Solidus surface: {ls_solidus_poly.GetNumberOfCells()} cells")

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

            # Solidus outline -- fixed #222222, solid line (see
            # MELT_BOUNDARY_COLOR above for why not the per-cut dark_color).
            solidus_rim = slice_at_x(ls_solidus_contour, x_cut)
            solidus_rim_poly = servermanager.Fetch(solidus_rim)
            solidus_segments = extract_polylines(solidus_rim_poly)
            log(f"  Solidus outline: {len(solidus_segments)} segments "
                f"(may be empty if no solid/mushy boundary straddles this x/window)")
            for seg in solidus_segments:
                add_line([x_disp_um] * len(seg), seg[:, 2] * 1e3 + z_off_mm, seg[:, 1] * 1e6 + y_off_um,
                         color=MELT_BOUNDARY_COLOR, linewidth=MELT_BOUNDARY_LINE_WIDTH)

            # Liquidus outline -- same fixed color, dashed to stay
            # distinguishable from the solid solidus line above.
            liquidus_rim = slice_at_x(ls_liquidus_contour, x_cut)
            liquidus_rim_poly = servermanager.Fetch(liquidus_rim)
            liquidus_segments = extract_polylines(liquidus_rim_poly)
            log(f"  Liquidus outline: {len(liquidus_segments)} segments "
                f"(may be empty if no liquid straddles this x/window)")
            for seg in liquidus_segments:
                add_line([x_disp_um] * len(seg), seg[:, 2] * 1e3 + z_off_mm, seg[:, 1] * 1e6 + y_off_um,
                         color=MELT_BOUNDARY_COLOR, linewidth=MELT_BOUNDARY_LINE_WIDTH,
                         linestyle=MELT_BOUNDARY_LIQUIDUS_DASH)

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
# --view=cutaway -- EXPERIMENTAL/PROTOTYPE, not part of the normal 4-view
# pipeline or _render_stacked_video.sh. Exactly render_lateral's own
# technique (gm_contour surface shell, colored by x_coord, camera looking
# down +x) -- just with the visible x-range restricted to
# [xmin, CUTAWAY_X_PLANE] instead of the full domain (user request,
# 2026-08-16: "the same thing as the side view but with only the surface
# that is x < -0.025"). Since render_lateral's own camera looks down the
# x-axis from large +x toward -x, near-centerline material (higher x)
# normally renders in front of, and can occlude, whatever sits farther
# away (lower x -- including the vapor-depression cavity's own walls);
# clipping that near-centerline slab away lets whatever surface remains at
# x < CUTAWAY_X_PLANE show through unobstructed, with no new rendering
# technique needed at all.
#
# CUTAWAY_X_PLANE = -0.025mm matches render_transverse's own "blue"
# cross-section (X_CROSS_SECTIONS[0] = -25e-6) -- same cut plane, different
# view. (An earlier version of this view used a much more elaborate
# solid-volume half-domain clip at x=-0.25mm with an oblique camera --
# abandoned after the first test render showed no melt reaches that far
# off-centerline, and after the user clarified they actually wanted this
# simpler render_lateral-based technique at the existing -0.025mm cut
# instead.)
#
# Deliberately independent of render_transverse/render_lateral: separate
# constants, separate function, no shared state -- neither is touched by
# this addition.
# ═════════════════════════════════════════════════════════════════════════
def render_cutaway(foam_file, time_value, output_png, output_pvsm, y_min=None, y_max=None, x_plane=None,
                    grid=False, top_crop_frac=0.10, show_velocity=False, supersample=1, dotted_solidus=False,
                    highlights=None, show_rays=False, show_section_velocity=False):
    FIELD_COLOR = 'x_coord'
    # x_plane: optional override of the cut plane itself (default matches
    # render_transverse's existing "blue" cross-section, X_CROSS_SECTIONS[0]
    # = -25um) -- same override pattern as y_min/y_max below, added for the
    # same cutaway-collage caller which now wants a second set of frames at
    # -30um alongside the original -25um set (user request, 2026-08-18).
    CUTAWAY_X_PLANE = -0.025e-3 if x_plane is None else x_plane  # m
    # y_min/y_max: optional per-call override of the shared Y_DEPTH_MIN/MAX
    # crop, for callers (e.g. the cutaway-collage script) that need a wider
    # window than the standard cutaway/cutaway3panel pipeline -- same
    # "view-specific override, don't unify with the shared constant" idea
    # as render_xray's own XRAY_Y_DEPTH_MAX, just parameterized instead of
    # hardcoded since only this one new caller wants it (user request,
    # 2026-08-18). None (the default) preserves existing behavior exactly.
    y_min = Y_DEPTH_MIN if y_min is None else y_min
    y_max = Y_DEPTH_MAX if y_max is None else y_max

    reader = OpenFOAMReader(FileName=foam_file)
    # 'U' only when needed -- extra field to read/interpolate for every
    # other caller of this function otherwise (user request, 2026-08-21:
    # prototype tangential-velocity arrows on the liquid surface).
    reader.CellArrays = [FIELD_GM, FIELD_LS, 'U'] if (show_velocity or show_section_velocity) else [FIELD_GM, FIELD_LS]
    reader.Createcelltopointfiltereddata = 1
    reader.UpdatePipeline(time=time_value)
    log("reader loaded")

    merged = MergeBlocks(Input=reader)
    merged.UpdatePipeline(time=time_value)
    log("blocks merged")

    # Cut away the far side of the domain (x < -150um) before any further
    # processing -- same idea as the CUTAWAY_X_PLANE cut on the near side,
    # just at the opposite extreme (user request, 2026-08-21). Applied this
    # early (on `merged` itself, reassigned) so every downstream
    # xmin-derived quantity (feature's own crop, camera framing, the
    # x-color transfer function's blue endpoint) picks up the new,
    # narrower domain automatically instead of needing a second adjustment.
    CUTAWAY_X_FAR_LIMIT = -0.150e-3  # m
    far_clip = Clip(Input=merged)
    far_clip.ClipType = 'Plane'
    far_clip.ClipType.Origin = [CUTAWAY_X_FAR_LIMIT, 0.0, 0.0]
    far_clip.ClipType.Normal = [1.0, 0.0, 0.0]
    far_clip.Invert = 0  # keep x >= CUTAWAY_X_FAR_LIMIT -- verified against
                          # the logged bounds below, not just assumed
    far_clip.UpdatePipeline(time=time_value)
    merged = far_clip

    bounds = merged.GetDataInformation().GetBounds()
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    log(f"Domain bounds (after x<{CUTAWAY_X_FAR_LIMIT * 1e3:.3f}mm far-side cutaway): "
        f"x=[{xmin},{xmax}] y=[{ymin},{ymax}] z=[{zmin},{zmax}]")

    laser_table = _load_laser_time_vs_position(os.path.dirname(foam_file))
    z_window_min, z_window_max = Z_VIEW_MIN, Z_VIEW_MAX
    log(f"Fixed crop window: x<{CUTAWAY_X_PLANE * 1e3:.3f}mm "
        f"y=[{y_min * 1e3:.3f},{y_max * 1e3:.3f}]mm "
        f"z=[{z_window_min * 1e3:.3f},{z_window_max * 1e3:.3f}]mm")

    laser_z = _laser_z_at(laser_table, time_value)
    log(f"Laser z={laser_z * 1e3:.3f}mm at t={time_value}")

    gm_contour = Contour(Input=merged)
    gm_contour.ContourBy = ['POINTS', FIELD_GM]
    gm_contour.Isosurfaces = [ISO_THRESHOLD]
    gm_contour.UpdatePipeline(time=time_value)
    gm_poly = servermanager.Fetch(gm_contour)
    log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

    # Mushy-zone surface, restricted to inside the metal first -- same
    # Clip-then-Contour pattern render_lateral uses (see its own comment
    # for why Clip, not Threshold).
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

    ls_liquidus_contour = Contour(Input=metal_only)
    ls_liquidus_contour.ContourBy = ['POINTS', FIELD_LS]
    ls_liquidus_contour.Isosurfaces = [LS_TLIQUIDUS]
    ls_liquidus_contour.UpdatePipeline(time=time_value)
    ls_liquidus_poly = servermanager.Fetch(ls_liquidus_contour)
    log(f"Mushy-zone (liquidus) surface: {ls_liquidus_poly.GetNumberOfCells()} cells")

    # Melt-boundary outline at the cut plane itself, not x=0 -- x=0 sits
    # well outside the kept x<CUTAWAY_X_PLANE region now, so an x=0
    # outline would mark a location that isn't shown at all.
    #
    # Slice-then-CONTOUR (not contour-then-slice, which is what this used
    # to do): contouring FIELD_GM/FIELD_LS directly on the 2D cut-plane
    # slice gives the boundary curve in one well-conditioned step, instead
    # of first building the full 3D marching-cubes surface and then
    # cutting *that* tessellated mesh with a plane (vtkCutter through a 3D
    # mesh) -- near-tangent intersections through sliver triangles in that
    # 3D mesh were fragmenting the outline into small disconnected dots/
    # gaps (bug caught by user, 2026-08-21; the gaps were real, pre-dating
    # the later LS_LINE_WIDTH*2 change, just less visible at the old
    # thinner width). Same technique the phase-fill code below (and
    # render_transverse) already uses for this same cut plane, just not
    # applied to the outline curves until now.
    outline_field_slice = Slice(Input=merged)
    outline_field_slice.SliceType = 'Plane'
    outline_field_slice.SliceType.Origin = [CUTAWAY_X_PLANE, (y_min + y_max) / 2.0, (z_window_min + z_window_max) / 2.0]
    outline_field_slice.SliceType.Normal = [1.0, 0.0, 0.0]
    outline_field_slice.UpdatePipeline(time=time_value)

    outline_metal_clip = Clip(Input=outline_field_slice)
    outline_metal_clip.ClipType = None
    outline_metal_clip.Scalars = ['POINTS', FIELD_GM]
    outline_metal_clip.Value = ISO_THRESHOLD
    outline_metal_clip.Invert = 0
    outline_metal_clip.UpdatePipeline(time=time_value)

    gm_curve = Contour(Input=outline_field_slice)
    gm_curve.ContourBy = ['POINTS', FIELD_GM]
    gm_curve.Isosurfaces = [ISO_THRESHOLD]
    gm_curve.UpdatePipeline(time=time_value)

    ls_curve = Contour(Input=outline_metal_clip)
    ls_curve.ContourBy = ['POINTS', FIELD_LS]
    ls_curve.Isosurfaces = [LS_TSOLIDUS]
    ls_curve.UpdatePipeline(time=time_value)

    ls_liquidus_curve = Contour(Input=outline_metal_clip)
    ls_liquidus_curve.ContourBy = ['POINTS', FIELD_LS]
    ls_liquidus_curve.Isosurfaces = [LS_TLIQUIDUS]
    ls_liquidus_curve.UpdatePipeline(time=time_value)

    def _weld(curve):
        """vtkCleanPolyData point-merging: Contour builds each line segment
        independently per mesh cell, so two segments that are supposed to
        share an endpoint often reference numerically-close-but-distinct
        point instances instead of the same one. Thick Wireframe rendering
        then draws each segment as its own independently-capped stroke,
        leaving a visible sub-pixel gap at every such joint (bug caught by
        user, 2026-08-21 -- persisted even after the slice-then-contour fix
        above, since that fixed *fragmentation* from cutting a 3D mesh, not
        this separate *unwelded-vertex* rendering issue). AbsoluteTolerance
        1e-9m is far below the ~5um mesh resolution, so this only welds
        points that were meant to be identical, not genuinely distinct
        nearby mesh points."""
        c = Clean(Input=curve)
        c.ToleranceIsAbsolute = 1
        c.AbsoluteTolerance = 1e-9
        c.UpdatePipeline(time=time_value)
        return c

    gm_curve = _weld(gm_curve)
    ls_curve = _weld(ls_curve)
    ls_liquidus_curve = _weld(ls_liquidus_curve)

    def _clip_curve_to_window(curve):
        """Box-clip an already-2D cut-plane curve down to the visible
        y/z window -- same bounds the old _outline_at_cut used, just no
        slicing step here since curve is already flat at CUTAWAY_X_PLANE."""
        c = Clip(Input=curve)
        c.ClipType = 'Box'
        c.ClipType.Position = [CUTAWAY_X_PLANE - 1e-6, y_min, z_window_min]
        c.ClipType.Length = [2e-6, y_max - y_min, z_window_max - z_window_min]
        c.Invert = 1
        c.UpdatePipeline(time=time_value)
        return c

    ls_slice_clip = _clip_curve_to_window(ls_curve)
    ls_slice_poly = servermanager.Fetch(ls_slice_clip)
    log(f"Solidus outline at x={CUTAWAY_X_PLANE * 1e3:.3f}mm: {ls_slice_poly.GetNumberOfCells()} cells "
        f"(may be empty if no melt currently straddles this cut plane)")

    ls_liquidus_slice_clip = _clip_curve_to_window(ls_liquidus_curve)
    ls_liquidus_slice_poly = servermanager.Fetch(ls_liquidus_slice_clip)
    log(f"Liquidus outline at x={CUTAWAY_X_PLANE * 1e3:.3f}mm: {ls_liquidus_slice_poly.GetNumberOfCells()} cells "
        f"(may be empty if no melt currently straddles this cut plane)")

    # Gas/metal boundary outline at the cut plane too (user request,
    # 2026-08-16).
    gm_slice_clip = _clip_curve_to_window(gm_curve)
    gm_slice_poly = servermanager.Fetch(gm_slice_clip)
    log(f"Gas/metal outline at x={CUTAWAY_X_PLANE * 1e3:.3f}mm: {gm_slice_poly.GetNumberOfCells()} cells "
        f"(may be empty if x<CUTAWAY_X_PLANE has no metal at all in this window)")

    # Spatial crop: y/z fixed (shared window), x restricted to
    # [xmin, CUTAWAY_X_PLANE] -- the one change from render_lateral's own
    # feature clip, which uses the full [xmin,xmax] x-range instead.
    feature = Clip(Input=gm_contour)
    feature.ClipType = 'Box'
    feature.ClipType.Position = [xmin, y_min, z_window_min]
    feature.ClipType.Length = [CUTAWAY_X_PLANE - xmin, y_max - y_min, z_window_max - z_window_min]
    feature.Invert = 1
    feature.UpdatePipeline(time=time_value)
    feature_poly = servermanager.Fetch(feature)
    log(f"Cropped feature: {feature_poly.GetNumberOfCells()} cells, "
        f"bounds={feature.GetDataInformation().GetBounds()}")

    xcolor = Calculator(Input=feature)
    xcolor.AttributeType = 'Point Data'
    xcolor.ResultArrayName = FIELD_COLOR
    xcolor.Function = 'coordsX*1e3'  # meters -> mm
    xcolor.UpdatePipeline(time=time_value)

    y_center = (y_min + y_max) / 2.0
    z_center = (z_window_min + z_window_max) / 2.0
    x_center = (xmin + xmax) / 2.0

    view = GetActiveViewOrCreate('RenderView')
    view.OrientationAxesVisibility = 0
    view.Background = [1, 1, 1]
    view.ViewSize = [max(1, round(VIEW_HEIGHT_PX * (z_window_max - z_window_min) / (y_max - y_min))), VIEW_HEIGHT_PX]
    view.ViewTime = time_value

    # highlights: a list of (z0,z1,y0,y1,color) sub-boxes, used below to
    # recolor the matching portions of the gm_outline curve (not the
    # surface fill -- tried that first, user request, 2026-08-21: "that's
    # not what I want", wanted the solid+liquid cross-section's outline
    # *line* colored instead) that fall inside a panel's own ROI box(es)
    # in the collage. Several panels can each carry their own highlight(s)
    # in their own color (user request, 2026-08-21: red for topleft/
    # topright, white for bottomleft, green+white for bottomright's two
    # boxes) -- see _tube_outline_multi below.
    highlights = highlights or []

    disp = Show(xcolor, view)
    disp.Representation = 'Surface'
    ColorBy(disp, ('POINTS', FIELD_COLOR))
    ctf = GetColorTransferFunction(FIELD_COLOR)
    # Back to x (not render_top's y -- tried that, user request 2026-08-17:
    # reverted) -- y duplicated exactly what the top panel already shows
    # (height), adding no new information; x is the one dimension an
    # orthographic side view otherwise flattens away entirely, so it's the
    # more informative choice specifically for this panel.
    #
    # Explicit transition bounds -- -135um to -35um (was -125/-25, user
    # request, 2026-08-22, to match this collage's actual --x-plane-um=-35
    # cut -- see below). Originally -125/-25 (user request, 2026-08-17 --
    # prior attempts: +-0.12mm fit-to-visible-range, -0.1/+0.05mm,
    # -0.1/0.0mm) so TRANSITION_RED_MM landed exactly on the *default*
    # CUTAWAY_X_PLANE (-25um); the collage script overrides that default to
    # -35um though, so -25um was actually past the real cut plane and the
    # visible edge never quite reached full red. -35um now lands exactly on
    # it again, so the cut plane itself -- the reddest point actually
    # reached -- is fully saturated red, using the whole blue->red range
    # across the visible surface instead of stopping short of it. x_min_mm
    # (the true domain edge) stays the fully-saturated flat blue endpoint;
    # the outer red endpoint is a dummy value just past TRANSITION_RED_MM,
    # purely so RGBPoints' values stay strictly increasing -- real data
    # never goes past it.
    x_min_mm = xmin * 1e3
    TRANSITION_BLUE_MM = -0.135
    TRANSITION_RED_MM = -0.035
    ctf.RGBPoints = [
        x_min_mm,                 0.0, 0.0, 1.0,
        TRANSITION_BLUE_MM,       0.0, 0.0, 1.0,
        TRANSITION_RED_MM,        1.0, 0.0, 0.0,
        TRANSITION_RED_MM + 0.01, 1.0, 0.0, 0.0,
    ]


    # ── Prototype: tangential-velocity arrows on the liquid surface ──────
    # (user request, 2026-08-21 -- "let's see how it looks"). Liquid-only
    # (T>=liquidus) subset of `feature`, tangential component of U (full U
    # minus its component along the local surface normal, so only in-plane
    # "surface flow" shows, not material moving toward/away from the
    # camera through the surface), sparsely sampled so arrows don't
    # blanket the whole melt pool. ScaleFactor is picked from this frame's
    # own actual |tangential U| range (not a fixed guess) so arrows land
    # at a sane visible length regardless of how fast this particular
    # frame's flow is.
    # Shared arrow style constants -- both this surface-velocity block and
    # the cross-section-velocity block further below (user request,
    # 2026-08-22: "same size and format as the arrows ... shown at the
    # surface") build their Glyphs from these same three values, so the two
    # arrow families are guaranteed identical in length/density rather than
    # just visually similar.
    VELOCITY_STRIDE = round(88 / 1.5)  # "1.5 times as dense as the current
                                    # state" (user feedback, 2026-08-21) --
                                    # density scales as 1/stride, so stride
                                    # shrinks from 88
    VELOCITY_OFFSET = 5             # "different seed" -- MaskPoints has
                                    # no RandomMode in this ParaView
                                    # version, so a nonzero Offset (which
                                    # points a given OnRatio stride
                                    # starts counting from) is the
                                    # deterministic, reproducible
                                    # equivalent: same density, different
                                    # subset of points (user feedback,
                                    # 2026-08-21)
    TARGET_ARROW_LEN_M = 90e-6 * 0.25  # "reduce sizes to 25%" (user
                                    # feedback, 2026-08-21), was 90e-6

    if show_velocity:
        liquid_only = Clip(Input=feature)
        liquid_only.ClipType = None
        liquid_only.Scalars = ['POINTS', FIELD_LS]
        liquid_only.Value = LS_TLIQUIDUS
        liquid_only.Invert = 0
        liquid_only.UpdatePipeline(time=time_value)

        # Clip outputs vtkUnstructuredGrid; GenerateSurfaceNormals requires
        # vtkPolyData, so re-extract the outer surface first.
        liquid_surface = ExtractSurface(Input=liquid_only)
        liquid_surface.UpdatePipeline(time=time_value)

        normals = GenerateSurfaceNormals(Input=liquid_surface)
        normals.UpdatePipeline(time=time_value)

        tang = Calculator(Input=normals)
        tang.AttributeType = 'Point Data'
        tang.ResultArrayName = 'TangentialU'
        tang.Function = 'U - (U.Normals)*Normals'
        tang.UpdatePipeline(time=time_value)

        tang_poly = servermanager.Fetch(tang)
        n_liquid_pts = tang_poly.GetNumberOfPoints()
        log(f"Liquid surface for velocity arrows: {n_liquid_pts} points")

        # Exactly LIQUID_FILL_COLOR (defined later in this function, for the
        # cut-plane phase fill) -- user request, 2026-08-21: "the same color
        # that we have for liquid at the cross-section". Duplicated here
        # (not referenced) since LIQUID_FILL_COLOR isn't defined yet at this
        # point in the function; keep the two literals in sync if either changes.
        VELOCITY_COLOR = [0.35, 0.35, 0.35]
        if n_liquid_pts > 0:
            tang_arr = tang_poly.GetPointData().GetArray('TangentialU')
            mags = np.linalg.norm(vtk_to_numpy(tang_arr), axis=1) if tang_arr is not None else np.array([0.0])
            max_mag = float(mags.max()) if mags.size else 0.0
            log(f"Tangential U: max |v|={max_mag:.4g} m/s")
        else:
            max_mag = 0.0

        if max_mag > 0:
            mask = MaskPoints(Input=tang)
            mask.OnRatio = VELOCITY_STRIDE
            mask.Offset = VELOCITY_OFFSET
            mask.GenerateVertices = 1
            mask.UpdatePipeline(time=time_value)

            # '2D Glyph'/Arrow -- flat line-drawn arrow, not the shaded 3D
            # cone+cylinder 'Arrow' source (user request, 2026-08-21: "not
            # a 3d object but simple line arrows").
            glyph = Glyph(Input=mask, GlyphType='2D Glyph')
            glyph.GlyphType.GlyphType = 'Arrow'
            glyph.OrientationArray = ['POINTS', 'TangentialU']
            # Uniform arrow length, not scaled by |tangential U| -- direction
            # only, no magnitude-by-size (user request, 2026-08-21).
            glyph.ScaleArray = ['POINTS', 'No scale array']
            glyph.ScaleFactor = TARGET_ARROW_LEN_M
            glyph.GlyphMode = 'All Points'
            glyph.UpdatePipeline(time=time_value)

            glyph_disp = Show(glyph, view)
            glyph_disp.Representation = 'Surface'
            glyph_disp.ColorArrayName = ['POINTS', '']
            glyph_disp.AmbientColor = VELOCITY_COLOR
            glyph_disp.DiffuseColor = VELOCITY_COLOR
            glyph_disp.LineWidth = LS_LINE_WIDTH
            log(f"Velocity glyphs: uniform ScaleFactor={glyph.ScaleFactor:.4g} "
                f"({TARGET_ARROW_LEN_M * 1e6:.0f}um every arrow), stride={VELOCITY_STRIDE}")
        else:
            log("No liquid surface / zero velocity this frame -- skipping velocity glyphs")

    # Still safely in front of the visible (x<CUTAWAY_X_PLANE) content --
    # doesn't need to change just because the visible range shrank.
    x_marker = xmax + 0.02 * (xmax - xmin)
    # Outline curves render OUTLINE_FRONT_NUDGE closer to the camera than
    # the fills below (larger x = closer, since the camera sits at large
    # +x looking toward -x) so they draw crisply on top instead of
    # z-fighting with the coincident fill surfaces.
    OUTLINE_FRONT_NUDGE = 5e-6  # m

    # Solid/mushy/liquid FILLED cross-section at the cut plane (user
    # request, 2026-08-16) -- reuses render_transverse's own field_slice ->
    # metal_clip -> clip_temp (3-way split by FIELD_LS) pipeline almost
    # verbatim (this file, render_transverse's per-cut loop), just Show()n
    # directly with a flat color per phase instead of being fetched into
    # matplotlib polygons: this view is already real ParaView Show()/
    # Render(), so that fits its own pipeline instead of bridging two
    # different rendering techniques. Translated to x_marker like the
    # outline curves, so it renders in front of the main (x-colored)
    # surface -- wherever a phase is absent (e.g. gas/void reaching this
    # far off-centerline), there's simply nothing drawn there, and the
    # main surface shows through from behind.
    # Neutral gray shades (user request, 2026-08-16 -- was red, was blue
    # before that) -- no existing render_transverse hue family to borrow
    # for gray, so these are hand-picked rather than duplicated from a
    # module/local constant like the earlier color versions were. Mushy
    # and liquid swapped from render_transverse's own light=solid/
    # dark=mushy/medium=liquid convention (user request, 2026-08-17: make
    # liquid darker than mushy) -- solid stays lightest, liquid is now
    # darkest, mushy is the medium shade. Now that the outline curves
    # above are pure black (LS_COLOR), all 3 shades sit comfortably
    # lighter than the lines so the crisp black boundaries still read
    # clearly on top of the fills.
    SOLID_FILL_COLOR = [0.75, 0.75, 0.75]
    MUSHY_FILL_COLOR = [0.55, 0.55, 0.55]
    LIQUID_FILL_COLOR = [0.35, 0.35, 0.35]

    field_slice = Slice(Input=merged)
    field_slice.SliceType = 'Plane'
    field_slice.SliceType.Origin = [CUTAWAY_X_PLANE, y_center, z_center]
    field_slice.SliceType.Normal = [1.0, 0.0, 0.0]
    field_slice.UpdatePipeline(time=time_value)

    metal_clip = Clip(Input=field_slice)
    metal_clip.ClipType = None
    metal_clip.Scalars = ['POINTS', FIELD_GM]
    metal_clip.Value = ISO_THRESHOLD
    metal_clip.Invert = 0
    metal_clip.UpdatePipeline(time=time_value)

    def _clip_temp(input_poly, value, invert):
        """One scalar Clip on FIELD_LS -- Invert=0 keeps T>=value, Invert=1
        keeps T<value (same convention as the gas/metal clip above)."""
        c = Clip(Input=input_poly)
        c.ClipType = None
        c.Scalars = ['POINTS', FIELD_LS]
        c.Value = value
        c.Invert = invert
        c.UpdatePipeline(time=time_value)
        return c

    def _fill_at_cut(phase_result):
        """Box-clip to the exact cut plane/y-z window (same bounds
        _outline_at_cut uses), then translate to x_marker."""
        c = Clip(Input=phase_result)
        c.ClipType = 'Box'
        c.ClipType.Position = [CUTAWAY_X_PLANE - 1e-6, y_min, z_window_min]
        c.ClipType.Length = [2e-6, y_max - y_min, z_window_max - z_window_min]
        c.Invert = 1
        c.UpdatePipeline(time=time_value)
        t = Transform(Input=c)
        t.Transform = 'Transform'
        t.Transform.Translate = [x_marker, 0.0, 0.0]
        t.UpdatePipeline(time=time_value)
        return t

    # solid: T<solidus. mushy: solidus<=T<liquidus (2 chained clips).
    # liquid: T>=liquidus. Same 3-way split as render_transverse.
    for phase_name, phase_color, phase_result in (
        ('solid', SOLID_FILL_COLOR, _clip_temp(metal_clip, LS_TSOLIDUS, 1)),
        ('mushy', MUSHY_FILL_COLOR, _clip_temp(_clip_temp(metal_clip, LS_TSOLIDUS, 0), LS_TLIQUIDUS, 1)),
        ('liquid', LIQUID_FILL_COLOR, _clip_temp(metal_clip, LS_TLIQUIDUS, 0)),
    ):
        phase_final = _fill_at_cut(phase_result)
        phase_poly = servermanager.Fetch(phase_final)
        log(f"{phase_name.capitalize()} fill at x={CUTAWAY_X_PLANE * 1e3:.3f}mm: "
            f"{phase_poly.GetNumberOfCells()} cells "
            f"(may be empty if no {phase_name} metal straddles this cut plane)")
        phase_disp = Show(phase_final, view)
        phase_disp.Representation = 'Surface'
        phase_disp.ColorArrayName = ['POINTS', '']
        phase_disp.AmbientColor = phase_color
        phase_disp.DiffuseColor = phase_color

    # Cross-section velocity arrows: same idea as the surface-velocity
    # block above, but flattened onto the flat cut-face slice instead of
    # the exterior 3D surface, covering the combined mushy+liquid region
    # there (user request, 2026-08-22: "white arrow in liquid + mushy
    # region at the cross-section ... same size and format as the arrows
    # ... shown at the surface"). The cut face's own "surface normal" is
    # simply the x-axis (it's flat), so there's no need for
    # GenerateSurfaceNormals/a per-point Normals array here -- ParaView's
    # Calculator has a constant unit-vector iHat built in, so
    # 'U - (U.iHat)*iHat' zeroes U's out-of-plane (x) component the same
    # way the surface block's 'U - (U.Normals)*Normals' does with its own
    # (per-point) normal. Built from metal_clip/_clip_temp/_fill_at_cut,
    # already defined above for the phase fill, with one extra
    # OUTLINE_FRONT_NUDGE step forward so the flat arrow glyphs don't
    # z-fight with the flat fill quad they sit on top of (same reasoning as
    # the outline curves' own x_marker+OUTLINE_FRONT_NUDGE).
    if show_section_velocity:
        liquid_mushy_slice = _clip_temp(metal_clip, LS_TSOLIDUS, 0)
        liquid_mushy_final = _fill_at_cut(liquid_mushy_slice)
        section_nudge = Transform(Input=liquid_mushy_final)
        section_nudge.Transform = 'Transform'
        section_nudge.Transform.Translate = [OUTLINE_FRONT_NUDGE, 0.0, 0.0]
        section_nudge.UpdatePipeline(time=time_value)

        planar = Calculator(Input=section_nudge)
        planar.AttributeType = 'Point Data'
        planar.ResultArrayName = 'PlanarU'
        planar.Function = 'U - (U.iHat)*iHat'
        planar.UpdatePipeline(time=time_value)

        planar_poly = servermanager.Fetch(planar)
        n_section_pts = planar_poly.GetNumberOfPoints()
        log(f"Liquid+mushy cut-face region for section-velocity arrows: {n_section_pts} points")

        SECTION_VELOCITY_COLOR = [1.0, 1.0, 1.0]  # white -- user request, 2026-08-22
        if n_section_pts > 0:
            planar_arr = planar_poly.GetPointData().GetArray('PlanarU')
            section_mags = np.linalg.norm(vtk_to_numpy(planar_arr), axis=1) if planar_arr is not None else np.array([0.0])
            section_max_mag = float(section_mags.max()) if section_mags.size else 0.0
            log(f"Planar (cut-face) U: max |v|={section_max_mag:.4g} m/s")
        else:
            section_max_mag = 0.0

        if section_max_mag > 0:
            section_mask = MaskPoints(Input=planar)
            section_mask.OnRatio = VELOCITY_STRIDE
            section_mask.Offset = VELOCITY_OFFSET
            section_mask.GenerateVertices = 1
            section_mask.UpdatePipeline(time=time_value)

            section_glyph = Glyph(Input=section_mask, GlyphType='2D Glyph')
            section_glyph.GlyphType.GlyphType = 'Arrow'
            section_glyph.OrientationArray = ['POINTS', 'PlanarU']
            section_glyph.ScaleArray = ['POINTS', 'No scale array']
            section_glyph.ScaleFactor = TARGET_ARROW_LEN_M
            section_glyph.GlyphMode = 'All Points'
            section_glyph.UpdatePipeline(time=time_value)

            section_glyph_disp = Show(section_glyph, view)
            section_glyph_disp.Representation = 'Surface'
            section_glyph_disp.ColorArrayName = ['POINTS', '']
            section_glyph_disp.AmbientColor = SECTION_VELOCITY_COLOR
            section_glyph_disp.DiffuseColor = SECTION_VELOCITY_COLOR
            section_glyph_disp.LineWidth = LS_LINE_WIDTH
            log(f"Section-velocity glyphs: uniform ScaleFactor={section_glyph.ScaleFactor:.4g} "
                f"({TARGET_ARROW_LEN_M * 1e6:.0f}um every arrow), stride={VELOCITY_STRIDE}")
        else:
            log("No liquid/mushy cut-face region / zero velocity this frame -- skipping section-velocity glyphs")

    # 2x the shared LS_LINE_WIDTH, local to render_cutaway's own 3
    # cut-plane outline curves only (not render_lateral's use of the same
    # shared constant) -- user request, 2026-08-21: "the lines that are
    # made at the cross-section (surface, and solidus and liquidus contour
    # lines) also need to get twice thick".
    CUTAWAY_OUTLINE_LINE_WIDTH = LS_LINE_WIDTH * 2 * 0.25  # 25% of the
                                        # previous value (user feedback,
                                        # 2026-08-21: "too thick")
    DOTTED_SOLIDUS_COLOR = [0.1, 0.35, 0.95]  # blue, collage_v3 prototype

    # Render these 3 outline curves as actual 3D tubes, not flat Wireframe
    # lines (user request, 2026-08-21: "cut-offs" -- root cause turned out
    # to be plain OpenGL line rendering's lack of joins between adjacent
    # segments: Contour emits one independent 2-point line cell per mesh
    # edge crossing, and at a sharp bend each segment's flat end-cap
    # leaves a wedge-shaped gap at the joint -- worse at larger LineWidth,
    # which is why doubling it made this newly visible. Point-welding
    # (vtkCleanPolyData, tried first) does NOT fix this: it only merges
    # coincident point *positions*, it doesn't turn independent line cells
    # into connected polylines, so the flat per-segment caps were still
    # there. A Tube filter sidesteps the whole issue -- real cylindrical
    # geometry with its own end caps, no dependence on neighboring cells
    # sharing a rendering primitive. World-space Radius is computed from
    # the same pixel width the old LineWidth used, via this view's own
    # world-units-per-pixel ratio, so it still looks like "an N-pixel-wide
    # line" on screen regardless of supersample.
    world_per_px = (y_max - y_min) * FRAME_MARGIN / VIEW_HEIGHT_PX
    outline_tube_radius = (CUTAWAY_OUTLINE_LINE_WIDTH * world_per_px) / 2.0

    def _stripped_producer(curve):
        # Clip (used upstream in ls_slice_clip/gm_slice_clip's box-crop)
        # outputs vtkUnstructuredGrid; Tube/Glyph require vtkPolyData.
        surf = ExtractSurface(Input=curve)
        surf.UpdatePipeline(time=time_value)
        # Tube tried directly on `surf` first -- came out *worse*, a
        # jagged "shark-tooth" mess (2026-08-21): Contour emits one
        # independent 2-point line cell per mesh-edge crossing, never
        # merged into continuous polylines, so Tube built a separate tiny
        # capped cylinder per segment instead of one continuous tube.
        # vtkStripper (not exposed as a paraview.simple proxy in this
        # ParaView build, so called directly via raw vtk + fetch/
        # TrivialProducer) joins those into maximal-length polylines
        # first -- JoinContiguousSegmentsOn because our segments are only
        # coincident-point-adjacent (via the earlier _weld/Clean pass),
        # not already cell-adjacent.
        poly = servermanager.Fetch(surf)
        stripper = vtk.vtkStripper()
        stripper.SetInputData(poly)
        stripper.JoinContiguousSegmentsOn()
        stripper.Update()
        producer = servermanager.sources.TrivialProducer()
        producer.GetClientSideObject().SetOutput(stripper.GetOutput())
        producer.UpdatePipeline()
        return producer

    def _tube_outline(curve, color=LS_COLOR, radius=None):
        producer = _stripped_producer(curve)
        t = Tube(Input=producer)
        t.Radius = outline_tube_radius if radius is None else radius
        t.NumberofSides = 8
        t.Capping = 1
        t.UpdatePipeline(time=time_value)
        disp = Show(t, view)
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['POINTS', '']
        disp.AmbientColor = color
        disp.DiffuseColor = color
        # Flat, unlit color instead of VTK's default Phong shading (user
        # report, 2026-08-22: the highlighted outlines "share the same
        # style that top of them is lighter and bottom is darker as they
        # are 3D, I don't want that") -- Ambient=1/Diffuse=0 means the
        # rendered color comes entirely from AmbientColor (constant,
        # independent of the tube's surface normal/orientation to the
        # light) rather than DiffuseColor (which is what produces the
        # direction-dependent lighter/darker shading on curved geometry).
        disp.Ambient = 1.0
        disp.Diffuse = 0.0
        return disp

    def _tube_outline_multi(curve, base_color, highlight_specs):
        """Same as _tube_outline, but each (z0,z1,y0,y1,color) box in
        highlight_specs recolors (at 2x radius) the portion of the curve
        that falls inside it -- one panel can carry several highlight
        boxes, each its own color (user request, 2026-08-21: red/white/
        green across topleft/topright/bottomleft/bottomright). Processed
        one box at a time: each pass clips its own "inside" piece off of
        whatever curve remains, so multiple boxes never fight over the
        same segment. Whatever's left after all boxes gets base_color at
        the normal radius.

        curve's Transform.Translate is a RELATIVE offset (not an absolute
        reposition), so its actual x is CUTAWAY_X_PLANE + x_marker + NUDGE
        -- not just x_marker + NUDGE, which was the first (wrong) attempt
        here: the box's x-range missed the curve entirely, so `inside` was
        always empty."""
        curve_x = CUTAWAY_X_PLANE + x_marker + OUTLINE_FRONT_NUDGE
        remaining = curve
        for z0, z1, y0, y1, color in highlight_specs:
            hl_position = [curve_x - 1e-6, y0, z0]
            hl_length = [2e-6, y1 - y0, z1 - z0]

            inside = Clip(Input=remaining)
            inside.ClipType = 'Box'
            inside.ClipType.Position = hl_position
            inside.ClipType.Length = hl_length
            inside.Invert = 1  # keep inside the box
            inside.UpdatePipeline(time=time_value)

            outside = Clip(Input=remaining)
            outside.ClipType = 'Box'
            outside.ClipType.Position = hl_position
            outside.ClipType.Length = hl_length
            outside.Invert = 0  # keep outside the box
            outside.UpdatePipeline(time=time_value)

            inside_poly = servermanager.Fetch(inside)
            log(f"Highlight: y=[{y0 * 1e6:.0f},{y1 * 1e6:.0f}]um "
                f"z=[{z0 * 1e3:.3f},{z1 * 1e3:.3f}]mm color={color} -- "
                f"{inside_poly.GetNumberOfPoints()} pts highlighted "
                f"(0 means the box missed the curve -- check it against the panel's actual box coords)")
            # A 'black' (== base_color/LS_COLOR) box is an *exclude* --
            # reverting a stray fragment back to looking like normal,
            # unhighlighted outline -- not an actual highlight, so it stays
            # at the normal 1x radius instead of the 2x used for real
            # colored highlights (user report, 2026-08-22: excluded regions
            # were rendering visibly thicker than the surrounding outline).
            is_exclude = color == [0.0, 0.0, 0.0]
            _tube_outline(inside, color, radius=outline_tube_radius if is_exclude else outline_tube_radius * 2)
            remaining = outside
        _tube_outline(remaining, base_color)

    def _dotted_outline(curve, color, dot_stride=6, dot_radius_factor=2.0):
        """Small sphere glyphs spaced along the curve instead of a solid
        tube -- a Tube has no dash-pattern equivalent (it's real 3D
        geometry, not a 2D stroke), so "dotted" is built literally as a
        series of dots (prototype, user request, 2026-08-21: dotted blue
        solid-liquid interface for collage_v3)."""
        producer = _stripped_producer(curve)
        mask = MaskPoints(Input=producer)
        mask.OnRatio = dot_stride
        mask.GenerateVertices = 1
        mask.UpdatePipeline(time=time_value)
        glyph = Glyph(Input=mask, GlyphType='Sphere')
        glyph.GlyphType.Radius = outline_tube_radius * dot_radius_factor
        glyph.ScaleArray = ['POINTS', 'No scale array']
        glyph.ScaleFactor = 1.0
        glyph.GlyphMode = 'All Points'
        glyph.UpdatePipeline(time=time_value)
        disp = Show(glyph, view)
        disp.Representation = 'Surface'
        disp.ColorArrayName = ['POINTS', '']
        disp.AmbientColor = color
        disp.DiffuseColor = color
        return disp

    ls_outline = Transform(Input=ls_slice_clip)
    ls_outline.Transform = 'Transform'
    ls_outline.Transform.Translate = [x_marker + OUTLINE_FRONT_NUDGE, 0.0, 0.0]
    ls_outline.UpdatePipeline(time=time_value)
    if dotted_solidus:
        # "Solid-liquid interface" = the solidus curve (T=Tsolidus, where
        # solid metal ends) -- distinct from the liquidus curve just below
        # it, which stays a normal black tube.
        ls_outline_disp = _dotted_outline(ls_outline, DOTTED_SOLIDUS_COLOR)
    else:
        ls_outline_disp = _tube_outline(ls_outline)

    # Liquidus now also LS_COLOR (black), not LS_LIQUIDUS_COLOR (gray) --
    # user request, 2026-08-16: solidus and liquidus both black. The 3
    # outline curves (gas/metal, solidus, liquidus) stay visually
    # distinguishable by their different y-depths, not by color anymore.
    ls_liquidus_outline = Transform(Input=ls_liquidus_slice_clip)
    ls_liquidus_outline.Transform = 'Transform'
    ls_liquidus_outline.Transform.Translate = [x_marker + OUTLINE_FRONT_NUDGE, 0.0, 0.0]
    ls_liquidus_outline.UpdatePipeline(time=time_value)
    ls_liquidus_outline_disp = _tube_outline(ls_liquidus_outline)

    gm_outline = Transform(Input=gm_slice_clip)
    gm_outline.Transform = 'Transform'
    gm_outline.Transform.Translate = [x_marker + OUTLINE_FRONT_NUDGE, 0.0, 0.0]
    gm_outline.UpdatePipeline(time=time_value)
    if highlights:
        _tube_outline_multi(gm_outline, LS_COLOR, highlights)
    else:
        gm_outline_disp = _tube_outline(gm_outline)

    # Camera: identical to render_lateral -- looking down +x, up = -y.
    view.CameraParallelProjection = 1
    view.CameraViewUp = [0, -1, 0]
    view.CameraFocalPoint = [x_center, y_center, z_center]
    view.CameraPosition = [x_center + 2.0 * (xmax - xmin), y_center, z_center]
    view.CameraParallelScale = (y_max - y_min) / 2.0 * FRAME_MARGIN
    Render(view)

    # supersample: render at a multiple of view.ViewSize (ParaView scales
    # line widths/fonts to match, not just pixel count) rather than at
    # view.ViewSize itself -- fixes jagged/low-quality outline curves and
    # (new, 2026-08-21) velocity arrows once the collage crops in tight on
    # a sub-region; default 1 preserves every existing caller's output
    # pixel-for-pixel (user request, 2026-08-21: "resolution ... too low").
    render_res = [round(view.ViewSize[0] * supersample), round(view.ViewSize[1] * supersample)]
    SaveScreenshot(output_png, view, ImageResolution=render_res)
    log(f"Saved: {output_png} at {render_res} ({supersample}x)")

    # Flatly overwrite the bottom 10% of the *raw* render with solid gray
    # to mask a rendering artifact in the bottom-right corner (user report,
    # 2026-08-22) -- done first, before any other post-processing (top
    # crop/grid/rays/colorbar) touches the image, since those only add
    # overlays or trim the *top* and would otherwise leave the corner
    # artifact untouched underneath.
    #
    # Fill color is *sampled* from the image itself (just above the fill
    # line, near the left edge -- the artifact is bottom-*right* only, so
    # this stays clean) rather than using the literal SOLID_FILL_COLOR
    # constant: VTK's default Phong shading renders that Ambient/
    # DiffuseColor triple visibly lighter on-screen than its raw RGB
    # values, so a hardcoded fill left a mismatched band (caught by eye
    # once rendered, 2026-08-22) -- sampling guarantees a pixel-exact match
    # regardless of how the renderer actually shades that flat fill.
    img = mpimg.imread(output_png)
    h_raw, w_raw = img.shape[:2]
    fill_start_row = round(h_raw * 0.90)
    sample_row = max(0, fill_start_row - 5)
    sample_col = round(w_raw * 0.05)
    fill_color = img[sample_row, sample_col, :3].copy()
    img[fill_start_row:, :, :3] = fill_color
    if img.shape[2] == 4:
        img[fill_start_row:, :, 3] = 1.0
    mpimg.imsave(output_png, img)
    log(f"Filled bottom 10% ({h_raw - fill_start_row}px) with sampled gray {fill_color} to mask corner artifact")

    CROP_FRAC = top_crop_frac
    if CROP_FRAC:
        _clip_top_fraction(output_png, CROP_FRAC)

    if grid or show_rays:
        # SURFACE_Y matches render_top's/main()'s own local constant.
        # Re-derives the extent this frame was actually rendered at from
        # the same camera/crop values just used above (view.ViewSize,
        # view.CameraParallelScale, CROP_FRAC) rather than guessing --
        # see _add_cutaway_grid's own docstring for why that matters. Shared
        # between --grid and --rays since both need to place an overlay in
        # this same physical-coordinate extent (user request, 2026-08-22).
        SURFACE_Y = 0.2e-3
        w0, h0 = view.ViewSize
        cps_y = view.CameraParallelScale       # half-height, world units
        cps_z = cps_y * (w0 / h0)              # half-width, world units (same aspect ratio)
        cut = round(h0 * CROP_FRAC)
        # Pre-crop top/bottom world y (row 0 / row h0-1); cropping removes
        # rows [0,cut) from the top only, so only the top edge moves.
        y_top_raw = y_center - cps_y
        y_bottom_raw = y_center + cps_y
        y_top_raw_postcrop = y_top_raw + cut * (y_bottom_raw - y_top_raw) / (h0 - 1)
        z_left_mm, z_right_mm = (z_center - cps_z) * 1e3, (z_center + cps_z) * 1e3
        y_top_offset_um = (y_top_raw_postcrop - SURFACE_Y) * 1e6
        y_bottom_offset_um = (y_bottom_raw - SURFACE_Y) * 1e6

        if grid:
            _add_cutaway_grid(
                output_png,
                z_left_mm=z_left_mm, z_right_mm=z_right_mm,
                y_top_offset_um=y_top_offset_um, y_bottom_offset_um=y_bottom_offset_um,
                time_us=time_value * 1e6,
            )

        if show_rays:
            _overlay_rays_on_cutaway(
                output_png, os.path.dirname(foam_file), time_value, ymin,
                z_left_mm=z_left_mm, z_right_mm=z_right_mm,
                y_top_offset_um=y_top_offset_um, y_bottom_offset_um=y_bottom_offset_um,
                supersample=supersample,
            )

    # Colorbar's title+bar rendered as usual (top-right placement,
    # labels_above=False (trying it back at the default -- below the bar,
    # like every other caller -- user request, 2026-08-17, to compare
    # against the labels_above=True version tried earlier the same day).
    # skip_overlay=True means it's saved as a standalone file only, not
    # pasted onto output_png itself -- this colorbar is meant to visually
    # overlap into the *top* panel above once stacked ("that moves it out
    # of the lateral-cutaway view and slightly into the top view but that
    # is cool", user request, 2026-08-17), which this function has no way
    # to do (it only ever sees this one panel's own canvas). The stacking
    # script (_tmp_stack3.sh) does that overlay instead, once the panels
    # are already combined -- see skip_overlay's own docstring.
    #
    # 2 labels (endpoints only), not 3 -- with values this close together
    # ("-0.125"/"-0.075"/"-0.025" are 6 characters each), a 3rd/middle
    # label overlapped its neighbors at this bar width and was illegible.
    #
    # Separate CTF for the colorbar's own visual range (+-CBAR_MARGIN_MM
    # around the transition bounds), not `ctf` itself (whose range extends
    # all the way to the true domain edge x_min_mm for actual data
    # coloring) -- same reasoning/technique as render_top's own colorbar
    # fix, user request, 2026-08-17.
    #
    # Labeled in um, not mm (user request, 2026-08-17) -- same
    # UM_PER_MM-scaling technique as render_top's own colorbar; harmless
    # here for the same reason (cbar_ctf is never bound to any real
    # Show()/ColorBy).
    UM_PER_MM = 1000.0
    cbar_ctf = GetColorTransferFunction(FIELD_COLOR + '_cbar_display')
    cbar_ctf.RGBPoints = [
        (TRANSITION_BLUE_MM - CBAR_MARGIN_MM) * UM_PER_MM, 0.0, 0.0, 1.0,
        TRANSITION_BLUE_MM * UM_PER_MM,                    0.0, 0.0, 1.0,
        TRANSITION_RED_MM * UM_PER_MM,                     1.0, 0.0, 0.0,
        (TRANSITION_RED_MM + CBAR_MARGIN_MM) * UM_PER_MM,  1.0, 0.0, 0.0,
    ]
    _overlay_colorbar(cbar_ctf, 'x (μm)', output_png,
                       custom_labels=[TRANSITION_BLUE_MM * UM_PER_MM, TRANSITION_RED_MM * UM_PER_MM],
                       labels_above=False, skip_overlay=True, thickness_scale=1.5)

    if output_pvsm:
        SaveState(output_pvsm)
        log(f"Saved state: {output_pvsm}")


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
        description="Render one of five VDEP power-sweep post-processing views.")
    parser.add_argument('--view', required=True, choices=['top', 'lateral', 'xray', 'transverse', 'cutaway'])
    parser.add_argument('case_foam')
    parser.add_argument('time', type=float)
    parser.add_argument('output_png')
    parser.add_argument('output_pvsm', nargs='?', default=None)
    parser.add_argument('--y-min-um', type=float, default=None,
                         help="cutaway view only: override the y-crop lower bound "
                              "(um, relative to nominal surface; negative = above "
                              "surface). Defaults to the shared Y_DEPTH_MIN if omitted.")
    parser.add_argument('--y-max-um', type=float, default=None,
                         help="cutaway view only: override the y-crop upper bound "
                              "(um, relative to nominal surface). Defaults to the "
                              "shared Y_DEPTH_MAX if omitted.")
    parser.add_argument('--x-plane-um', type=float, default=None,
                         help="cutaway view only: override the cut plane (um). "
                              "Defaults to -25 (matches render_transverse's "
                              "X_CROSS_SECTIONS[0]) if omitted.")
    parser.add_argument('--grid', action='store_true',
                         help="cutaway view only: also save a *_grid.png copy with a "
                              "z(mm)/y-offset-from-surface(um) coordinate grid overlaid, "
                              "as an annotation-planning aid. Original output_png is untouched.")
    parser.add_argument('--top-crop-frac', type=float, default=0.10,
                         help="cutaway view only: fraction of the raw render's top rows to "
                              "discard (default 0.10, matching the standard cutaway/cutaway3panel "
                              "pipeline). Pass 0 to keep the full requested y-window uncropped -- "
                              "needed when the top of the y-range is meant to show as blank white "
                              "(e.g. the collage panels' y-min sitting above the domain's real "
                              "edge) rather than being silently trimmed away.")
    parser.add_argument('--velocity', action='store_true',
                         help="cutaway view only: overlay tangential-velocity arrows (green) on "
                              "the liquid portion of the surface -- prototype, 2026-08-21.")
    parser.add_argument('--supersample', type=int, default=1,
                         help="cutaway view only: render at this multiple of the normal pixel "
                              "resolution (ParaView scales line widths/fonts to match), for "
                              "sharper outline curves/arrows once cropped in tight for a collage. "
                              "Default 1 (no change).")
    parser.add_argument('--dotted-solidus', action='store_true',
                         help="cutaway view only: draw the solidus (solid-liquid interface) curve "
                              "as a dotted blue line of small spheres instead of a solid black tube "
                              "-- prototype, 2026-08-21, for collage_v3.")
    parser.add_argument('--highlight', action='append', default=[],
                         help="cutaway view only: recolor the gm_outline curve inside a (z,y) box. "
                              "Repeatable -- one panel can carry several, each its own color. "
                              "Format: z0,z1,y0,y1,color -- z in mm, y as an offset from the "
                              "nominal surface in um (same convention as --y-min-um/--y-max-um), "
                              "color one of red/white/green.")
    parser.add_argument('--rays', action='store_true',
                         help="cutaway view only: overlay laser ray-tracing segments (dark red, "
                              "0-25%% opacity scaled by each segment's power) from the solver's "
                              "VTKs/rays_laser0.vtk.series, if this case has one -- prototype, "
                              "2026-08-22, for collage_v3. No-ops (with a log line) if the case "
                              "has no ray VTK series.")
    parser.add_argument('--section-velocity', action='store_true',
                         help="cutaway view only: overlay flow-direction arrows (white) across the "
                              "combined mushy+liquid region of the flat cut-face slice, same size/"
                              "density/style as --velocity's surface arrows -- for collage_v4, "
                              "2026-08-22.")
    args = parser.parse_args()

    if args.view == 'top':
        render_top(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'lateral':
        render_lateral(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'transverse':
        render_transverse(args.case_foam, args.time, args.output_png, args.output_pvsm)
    elif args.view == 'cutaway':
        # SURFACE_Y matches render_top's own local constant -- nominal
        # flat-plate surface height, see topoSetDict's "y surface (0.2mm)"
        # and setFieldsDict's "metal: y=0.2mm to 0.5mm; 0.2mm gas above".
        # --y-min-um/--y-max-um are offsets from this surface (negative =
        # above surface/gas side, positive = below surface/metal side,
        # same sign convention as render_top's own y-offset coloring) --
        # NOT raw domain y, which is what these silently did before this
        # fix (bug caught by user, 2026-08-18: a requested +-300um window
        # rendered mostly blank because it was centered on raw y=0, which
        # is actually 200um *above* the true surface, not on the surface
        # itself).
        SURFACE_Y = 0.2e-3
        y_min = SURFACE_Y + args.y_min_um * 1e-6 if args.y_min_um is not None else None
        y_max = SURFACE_Y + args.y_max_um * 1e-6 if args.y_max_um is not None else None
        x_plane = args.x_plane_um * 1e-6 if args.x_plane_um is not None else None
        HIGHLIGHT_COLORS = {
            'red': [1.0, 0.639, 0.867],         # matches the collage script's PROTRUSION_LABEL_COLOR (#FFA3DD,
                                                 # a light pink -- user tried this directly, 2026-08-22, after the
                                                 # crimson #e6005a fixed the same low-contrast complaint)
            'white': [1.0, 1.0, 1.0],
            'green': [0.7176, 0.9608, 0.5569],  # matches the collage script's PORE0_LABEL_COLOR (#b7f58e)
            'blue': [0.5569, 0.8039, 0.9608],   # matches the collage script's CONTACT_LABEL_COLOR/BOX_COLOR (#8ecdf5)
            'pore2': [0.5838, 0.9608, 0.5569],  # matches the collage script's PORE1_LABEL_COLOR (#95f58e)
            'black': [0.0, 0.0, 0.0],           # explicit "un-highlight"/exclude; must be listed
                                                 # BEFORE the larger box it should carve out of
        }
        highlights = []
        for spec in args.highlight:
            z0_mm, z1_mm, y0_um, y1_um, color = spec.split(',')
            highlights.append((
                float(z0_mm) * 1e-3, float(z1_mm) * 1e-3,
                SURFACE_Y + float(y0_um) * 1e-6, SURFACE_Y + float(y1_um) * 1e-6,
                HIGHLIGHT_COLORS[color],
            ))
        render_cutaway(args.case_foam, args.time, args.output_png, args.output_pvsm,
                        y_min=y_min, y_max=y_max, x_plane=x_plane, grid=args.grid,
                        top_crop_frac=args.top_crop_frac, show_velocity=args.velocity,
                        supersample=args.supersample, dotted_solidus=args.dotted_solidus,
                        highlights=highlights, show_rays=args.rays,
                        show_section_velocity=args.section_velocity)
    elif args.view == 'xray':
        if args.output_pvsm:
            log("Note: --view=xray never had a ParaView state to save; ignoring output_pvsm arg.")
        render_xray(args.case_foam, args.time, args.output_png)


if __name__ == '__main__':
    main()
