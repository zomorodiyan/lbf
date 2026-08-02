# Transverse cross-section view of the melt pool -- normal ParaView
# screenshot variant, three panels at fixed distances behind the laser.
#
# Sibling to lateral_screenshot.py/top_screenshot.py: same technique
# (marching-cubes isosurface of alpha_smoothed, rendered directly by
# ParaView rather than ray-traced) -- but "transverse" here means
# perpendicular to the scan direction: looking straight down the z-axis
# (scan/travel axis) at the x-y plane (lateral width vs. depth), the classic
# melt-pool cross-section shot from AM literature. This is a different plane
# from lateral_screenshot.py's "lateral" view, which looks down x (through-
# thickness) at the y-z (depth vs. scan-length) plane -- effectively a
# *longitudinal* section, despite the name.
#
# A cross-section needs an actual cut, not a full-domain render: looking
# straight down z at the whole opaque melt-track surface would just show
# whichever surface point is nearest the camera (the current laser-adjacent
# frontier), occluding everything behind it -- not a cross-section at any
# particular position. So each panel clips the gas/metal isosurface at a
# target cut location and renders what's behind it -- same "clip then render
# the isosurface directly" technique as the other two scripts' fixed crop
# windows, just with the z-window now bounded in front by the cut instead of
# spanning the whole scan length.
#
# The clip's *back* edge is the full domain start (zmin), not a thin slab:
# with an opaque, orthographic, straight-down--z camera, a fragment only
# shows up in the final image if it's the nearest-to-camera (largest-z) hit
# at its (x,y) -- so extending the kept region further back than whatever's
# already nearest at each pixel can never change what's drawn there. What it
# *does* change is every pixel where the near-cut region has no track
# geometry at all (most of the frame, away from the melt track): instead of
# an empty gap (nothing kept that far back) or a knife-edge sliver of the
# flat plate (only a thin slab kept), the actual flat plate/powder bed
# extending all the way back renders as the solid, continuous part it really
# is. Hence "the part is fully there for z < cut" -- this is a strictly more
# complete picture of the same camera view, not a different one.
#
# Per-frame reference point: the laser's own current z position (from
# timeVsLaserPosition, piecewise-linear-interpolated -- same helper as
# lateral_xray.py), not a fixed z -- so the same set of "distance behind the
# laser" cuts stays meaningful at any timestep as the track grows, instead of
# needing a hand-picked fixed z per case/frame.
#
# Two nested surfaces per panel, both fully opaque:
#   1. Gas/metal surface (alpha_smoothed=0.5) -- the part's outer boundary.
#      Colored by dist_behind_laser (laser_z - coordsZ, in mm, computed via
#      Calculator) on a smooth REL_Z_COLD-REL_Z_HOT scale -- distance behind
#      the *laser's current position*, not absolute z and not depth-
#      relative-to-cut. Range: REL_Z_COLD is the first-listed panel's own
#      offset (also the furthest-from-laser one, since OFFSETS_BEHIND_LASER
#      is listed largest-first), REL_Z_HOT is the furthest panel offset plus
#      1mm -- so the whole set of panels sits inside one consistent,
#      comparable scale. (A two-tone solid/mushy gray split was tried here
#      instead and reverted -- this gradient scheme was fine as-is.)
#   2. Liquid/solid (mushy-zone) surface (T=TSolidus) -- entirely nested
#      inside surface 1 (mushy material only exists within metal, never
#      outside it). Flat light gray (LS_SURFACE_COLOR), not colored --
#      distinguishes it from surface 1 by color, not just position. Both
#      opaque, so surface 2 is only visible where it's the frontmost
#      (nearest-camera) geometry at a given pixel -- i.e. through gaps in
#      surface 1's own coverage, or at/near the rim -- not revealed via
#      transparency (tried earlier; reverted per instruction).
#
# Both surfaces additionally get an exact-plane rim outline at the cut
# (slice_at_cut, below) so their cross-sectional boundary is visible
# regardless of the above occlusion: black for the mushy-zone rim; the
# gas/metal rim gets one of GM_RIM_COLORS per panel/offset (not one shared
# color), so the same per-offset colors can later mark each cut's location
# on the top-view script (top_screenshot.py) too.
#
# Cross-section cap (cross_section_cap): unlike the rim outlines above (1D
# curves), this is an actual *filled* 2D cross-section -- Slice the metal
# volume itself (metal_only, not just its outer gm_contour shell) at z=z_cut,
# giving the true cross-sectional area of the part at that plane, then keep
# only the solid (T<TSolidus) portion of it and fill that in flat
# CAP_SOLID_COLOR (light gray). Sits exactly at the clip window's own front
# face, so it's the frontmost possible geometry -- the mushy portion of the
# cap is deliberately left unfilled, so the melt pool reads as a "hole" in
# the light-gray solid face, through which the colored/textured surfaces
# behind are visible -- the classic metallographic-cross-section look.
#
# (A "waterline" variant -- T=TSolidus contoured directly on the gm_contour
# surface mesh, embedded in the skin instead of a separate nested surface --
# was tried and reverted; didn't read well.)
#
# Mushy-zone field choice: T (temperature), not epsilon1. epsilon1 is
# *derived* from T inside the solver (TEqn.H computes it from an enthalpy
# correction, then clamps it to [0,1]), so it saturates to a hard 0 or 1
# almost everywhere and only genuinely varies across a very thin band --
# exactly the kind of near-step field that contours poorly. T is the
# primary, continuous field the energy equation actually solves for, with no
# clamping, so its solidus crossing is much better-conditioned for marching
# cubes. TLiquidus is also defined below if the fully-liquid boundary is
# wanted in addition to/instead of the solidus one currently used.
#
# Rim outlines are an exact geometric Slice (a plane intersection) at
# z=z_cut, not a rendered clip-and-show: this sidesteps the depth-tie
# problem entirely rather than just narrowing its window -- a Slice computes
# the exact 1D intersection curve mathematically, so there is no
# rendered-surface depth-buffer ambiguity left to arbitrate in the first
# place (an earlier fix -- clipping to a thin window right at the cut --
# only shrank how much of that ambiguity could occur; this removes the
# mechanism outright).
#
# Offsets chosen to land within the melt-affected zone: checked empirically
# (t=4e-4s, testrun64), the mushy zone only extends about the last ~0.5mm
# behind the laser, so OFFSETS_BEHIND_LASER stays within that range -- unlike
# an earlier version of this script (1.0/1.5/2.0mm), which was well past it
# and always read empty for the mushy-zone surface/rim. Still time/case-
# dependent, though: at other timesteps or cases the melt-affected distance
# may differ, in which case these offsets could again read empty (correctly
# reflecting no melt left there, not a bug) or need adjusting.
#
# Compositing: each of the 3 offsets is rendered to its own temporary PNG
# with the exact same per-panel pipeline (crop, color, camera -- bare
# render, no colorbar/text overlay, see below), then stacked side-by-side
# into the final
# output with matplotlib/numpy -- simpler and more robust than driving
# multiple ParaView RenderViews inside one Layout for a single screenshot.
# All panels share one fixed crop window (unlike the other two scripts, whose
# z-window spans the whole scan and so needs a growing-track trim step),
# so no post-hoc whitespace trim is needed here.
#
# Run via pvpython (needs paraview.simple to read the OpenFOAM case):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/transverse_screenshot.py /workspace/<case.foam> <time> <output.png> [<output.pvsm>]
#
from paraview.simple import *
from paraview import servermanager
import paraview.simple
import sys
import os
import re
import time
import tempfile
import shutil

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
FIELD_LS      = 'T'               # liquid/solid (mushy-zone) interface source
                                   # -- temperature, not epsilon1 (see header
                                   # for why: epsilon1 is a clamped/derived
                                   # field, T is the continuous primary one)
FIELD_COLOR   = 'dist_behind_laser'  # gas/metal surface coloring -- distance
                                   # behind the laser's current z position
                                   # (mm), computed below via Calculator; NOT
                                   # temperature, absolute z, or depth-
                                   # relative-to-cut (see header)
ISO_THRESHOLD = 0.5               # gas/metal isosurface value (alpha_smoothed).
LS_TSOLIDUS  = 840.0               # K -- mushy-zone bounds, per CLAUDE.md's
LS_TLIQUIDUS = 867.0               # K -- physical parameters (AlSi10Mg)
CAP_SOLID_COLOR = [0.8, 0.8, 0.8]   # flat light gray for the solid portion of
                                     # the filled cross-section cap (see header)
LS_SURFACE_COLOR = [0.75, 0.75, 0.75]  # flat lighter gray for the mushy-zone
                                        # surface itself (was 0.5,0.5,0.5) --
                                        # not colored
LS_COLOR = [0.0, 0.0, 0.0]         # flat color for the mushy-zone rim outline
                                    # -- bold black (pink tried and reverted,
                                    # didn't read well)
LS_LINE_WIDTH = 4.0                 # wireframe line width (px) for the outline
GM_RIM_COLORS = [                    # one flat color per offset/panel for the
    [1.0, 0.0, 0.0],                 # gas/metal rim outline -- distinct per
    [1.0, 0.5, 0.0],                 # cross-section (not one shared color)
    [1.0, 1.0, 0.0],                 # so the same colors can mark each cut's
]                                    # location on the top-view script too.
                                     # Indexed in OFFSETS_BEHIND_LASER order:
                                     # red, orange, yellow.
GM_RIM_LINE_WIDTH = 3.0              # wireframe line width (px) for the rim
OFFSETS_BEHIND_LASER = [0.7e-3, 0.6e-3, 0.5e-3]  # meters behind the laser --
                                                   # user-chosen (2026-08-01),
                                                   # red/orange/yellow = 0.7/
                                                   # 0.6/0.5mm. Sits right at
                                                   # the edge of, or just past,
                                                   # the ~0.5mm melt-affected
                                                   # zone found empirically in
                                                   # the header note above --
                                                   # so the mushy-zone surface/
                                                   # rim may read empty for
                                                   # some panels at some
                                                   # timesteps (correctly, not
                                                   # a bug -- see header).
                                                   # Largest offset (smallest,
                                                   # most-behind z) first --
                                                   # left-to-right panel order
                                                   # increases in z, matching
                                                   # lateral_screenshot.py's
                                                   # left-to-right z convention
REL_Z_COLD = OFFSETS_BEHIND_LASER[0] * 1e3        # mm -- color-scale bounds
                                                    # for dist_behind_laser:
                                                    # the first-listed panel's
                                                    # own offset (also the
                                                    # largest/furthest-from-
                                                    # laser one, since offsets
                                                    # are listed largest-first)
REL_Z_HOT = max(OFFSETS_BEHIND_LASER) * 1e3 + 1.0  # ...to the furthest
                                                     # panel's offset + 1mm
X_LATERAL_MIN = -0.2e-3           # x crop (screen-horizontal) -- matches
X_LATERAL_MAX = 0.2e-3             # top_screenshot.py's crop window
Y_DEPTH_MIN = 0.05e-3              # y crop (screen-vertical) -- matches
Y_DEPTH_MAX = 0.4e-3                # lateral_screenshot.py's crop window
VIEW_HEIGHT_PX = 500               # per-panel output height in px; width is
                                    # set to match the crop window's own
                                    # x/y aspect ratio, so there's no
                                    # letterboxing
FRAME_MARGIN = 1.3                 # small margin so content doesn't fill the
                                    # frame edge-to-edge. This value's own
                                    # original justification ("the colorbar/
                                    # text overlay has room") was stale --
                                    # this script's panels are a bare render
                                    # with no colorbar or text overlay at all
                                    # (see "No colorbar/text label overlay"
                                    # below) -- but left at 1.3 regardless:
                                    # unlike lateral_screenshot.py/
                                    # top_screenshot.py, this script has no
                                    # post-render crop-to-content step to
                                    # begin with (see header: "no post-hoc
                                    # whitespace trim is needed here"), so
                                    # tightening this wouldn't shrink the
                                    # output at all, just change the zoom
                                    # level within the same fixed canvas
PANEL_GAP_PX = 20                  # white gap between panels in the composite
# ─────────────────────────────────────────────────────────────────

if len(sys.argv) not in (4, 5):
    print("Usage: pvpython transverse_screenshot.py <case.foam> <time> <output.png> [<output.pvsm>]")
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
laser_z = _laser_z_at(laser_table, time_value)
log(f"Laser z={laser_z*1e3:.3f}mm at t={time_value}")

gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

# Mushy-zone (liquid/solid) surface, restricted to inside the metal first --
# same metal_only-then-contour pattern as lateral_xray.py's ls_contour, so
# T isn't contoured across the (physically unrelated there) gas/plume region.
#
# Clip by scalar (not Threshold): metal_only feeds both ls_contour (the
# mushy surface/rim, compared against the exact gm_contour surface) and
# cap_full below (the cross-section cap, whose outer/gas-facing edge should
# align with gm_rim). Threshold's blocky, cell-based boundary here was
# exactly the same "cell artifact" class of mismatch already fixed for the
# cap's solid/mushy split -- just one level further upstream, affecting both
# the gas-vs-mushy-surface and gas-vs-cap-solid boundaries at once. Clip with
# ClipType=None clips by Scalars/Value directly, interpolated like Contour,
# matching gm_contour's own precision. Invert=0 confirmed empirically to
# match Threshold's old [ISO_THRESHOLD, 1.0] selection (identical bounds).
metal_only = Clip(Input=merged)
metal_only.ClipType = None
metal_only.Scalars = ['POINTS', FIELD_GM]
metal_only.Value = ISO_THRESHOLD
metal_only.Invert = 0
metal_only.UpdatePipeline(time=time_value)

ls_contour = Contour(Input=metal_only)
ls_contour.ContourBy = ['POINTS', FIELD_LS]
ls_contour.Isosurfaces = [LS_TSOLIDUS]  # single line: the outer edge of the
                                          # melt-affected zone (anything past
                                          # solidus has been at least
                                          # partially liquid) -- the standard
                                          # "melt pool boundary" definition.
                                          # LS_TLIQUIDUS is still defined
                                          # above if the fully-liquid line is
                                          # wanted instead/also.
ls_contour.UpdatePipeline(time=time_value)
ls_poly = servermanager.Fetch(ls_contour)
log(f"Mushy-zone (liquid/solid) surface: {ls_poly.GetNumberOfCells()} cells")

x_center = (X_LATERAL_MIN + X_LATERAL_MAX) / 2.0
y_center = (Y_DEPTH_MIN + Y_DEPTH_MAX) / 2.0
view_w = max(1, round(VIEW_HEIGHT_PX * (X_LATERAL_MAX - X_LATERAL_MIN) / (Y_DEPTH_MAX - Y_DEPTH_MIN)))


def box_clip(input_poly, z_window_min, z_window_max):
    """x/y fixed crop window (see header), z in [z_window_min, z_window_max]."""
    c = Clip(Input=input_poly)
    c.ClipType = 'Box'
    c.ClipType.Position = [X_LATERAL_MIN, Y_DEPTH_MIN, z_window_min]
    c.ClipType.Length = [X_LATERAL_MAX - X_LATERAL_MIN, Y_DEPTH_MAX - Y_DEPTH_MIN, z_window_max - z_window_min]
    c.Invert = 1  # keep the inside of the box
    c.UpdatePipeline(time=time_value)
    return c


def slice_at_cut(input_poly, z_cut):
    """Exact plane intersection at z=z_cut (see header), then crop to the
    x/y window -- no z window needed since the slice is already flat at
    z_cut."""
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
            f"kept z=[{zmin*1e3:.3f},{z_cut*1e3:.3f}]mm (full part, see header)")

        # Spatial crop (see header): x/y fixed window, z from the full domain
        # start through the cut -- "the part, fully there" behind the cut,
        # not just a thin slab.
        feature = box_clip(gm_contour, zmin, z_cut)
        feature_poly = servermanager.Fetch(feature)
        log(f"  Cropped: {feature_poly.GetNumberOfCells()} cells, "
            f"bounds={feature.GetDataInformation().GetBounds()}")

        # Liquid/solid surface, same full-depth clip as the gas/metal one --
        # nested inside it (see header on why it's mostly hidden, opaque).
        ls_surface = box_clip(ls_contour, zmin, z_cut)
        ls_surface_poly = servermanager.Fetch(ls_surface)
        log(f"  Cropped mushy-zone surface: {ls_surface_poly.GetNumberOfCells()} cells "
            f"(see header -- may be empty this far behind the laser)")

        # Exact plane intersection at the cut, not a rendered clip -- see
        # header for why this replaces the earlier thin-box-clip approach.
        ls_feature = slice_at_cut(ls_contour, z_cut)
        ls_feature_poly = servermanager.Fetch(ls_feature)
        log(f"  Mushy-zone rim at cut: {ls_feature_poly.GetNumberOfCells()} cells "
            f"(see header -- may be empty this far behind the laser)")

        gm_rim = slice_at_cut(gm_contour, z_cut)
        gm_rim_poly = servermanager.Fetch(gm_rim)
        log(f"  Gas/metal rim outline at cut: {gm_rim_poly.GetNumberOfCells()} cells")

        # dist_behind_laser (mm) for coloring the gas/metal surface -- laser_z
        # minus each point's own z, not absolute z or depth-relative-to-cut
        # (see header).
        gm_zcolor = Calculator(Input=feature)
        gm_zcolor.AttributeType = 'Point Data'
        gm_zcolor.ResultArrayName = FIELD_COLOR
        gm_zcolor.Function = f'({laser_z}-coordsZ)*1e3'
        gm_zcolor.UpdatePipeline(time=time_value)

        # Filled cross-section cap (see header): slice the metal *volume*
        # (not just its outer shell) at z=z_cut, then keep only the solid
        # portion. Genuinely different object from gm_rim/ls_feature above --
        # those are 1D curves from slicing surfaces; this is a 2D filled area
        # from slicing the volume.
        cap_full = slice_at_cut(metal_only, z_cut)
        # Clip by scalar value (T), not Threshold: Threshold keeps/discards
        # whole cells, giving a blocky/staircase boundary at cell resolution
        # ("cell artifacts"). Clip with ClipType=None clips by the Scalars/
        # Value pair directly, cutting cells exactly where T crosses
        # LS_TSOLIDUS (interpolated, like Contour) instead of at cell
        # boundaries -- confirmed empirically to give a visibly smooth edge
        # vs. Threshold's jagged one.
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
        view.ViewTime = time_value  # see lateral_screenshot.py -- the view has
                                     # its own time state independent of the
                                     # per-filter UpdatePipeline(time=...) calls

        # Shown first (flat gray, not colored) -- distinguishes it from the
        # gas/metal surface by color, not just position. Both fully opaque,
        # so ordinarily the z-buffer alone decides which one wins at a given
        # pixel regardless of Show() order -- but gm is drawn *after* this,
        # so it wins any exact-tie/z-fighting case too (see header: gm is
        # always the true outer boundary, mush is always nested inside it,
        # so gm should never legitimately lose to it).
        ls_surface_disp = Show(ls_surface, view)
        ls_surface_disp.Representation = 'Surface'
        ls_surface_disp.ColorArrayName = ['POINTS', '']  # see ls_disp below
                                                           # for why not
                                                           # ColorBy(rep, None)
        ls_surface_disp.AmbientColor = LS_SURFACE_COLOR
        ls_surface_disp.DiffuseColor = LS_SURFACE_COLOR

        disp = Show(gm_zcolor, view)
        disp.Representation = 'Surface'
        ColorBy(disp, ('POINTS', FIELD_COLOR))

        ctf = GetColorTransferFunction(FIELD_COLOR)
        # Smooth blue->white gradient across the full REL_Z_COLD-REL_Z_HOT
        # range. White is the endpoint (at REL_Z_HOT), not a midpoint on the
        # way to red -- ParaView's default diverging interpolation between
        # two RGBPoints would otherwise wash out through white *then*
        # continue to whatever the second point's color is; setting that
        # second point itself to white stops the gradient there.
        ctf.RGBPoints = [
            REL_Z_COLD, 0.0, 0.0, 1.0,
            REL_Z_HOT,  1.0, 1.0, 1.0,
        ]

        # Cross-section cap's solid portion, drawn last (frontmost -- see
        # header: it sits exactly at the clip window's own front face, so
        # this only matters for exact ties, not real occlusion).
        cap_solid_disp = Show(cap_solid, view)
        cap_solid_disp.Representation = 'Surface'
        cap_solid_disp.ColorArrayName = ['POINTS', '']
        cap_solid_disp.AmbientColor = CAP_SOLID_COLOR
        cap_solid_disp.DiffuseColor = CAP_SOLID_COLOR

        # No colorbar/text label overlay -- bare render, matching
        # _probe_close_offset.png's look.
        ls_disp = Show(ls_feature, view)
        ls_disp.Representation = 'Wireframe'
        # Not ColorBy(ls_disp, None) -- when ls_feature is empty (no melt
        # this far behind the laser, common at early timesteps/far offsets),
        # Show() has no array to default to, and ColorBy(rep, None) then
        # crashes trying to look up an association for that "no array" state.
        # Setting ColorArrayName directly to an empty name is the same
        # "use solid color" end state without going through that lookup.
        ls_disp.ColorArrayName = ['POINTS', '']
        ls_disp.AmbientColor = LS_COLOR
        ls_disp.DiffuseColor = LS_COLOR
        ls_disp.LineWidth = LS_LINE_WIDTH

        gm_rim_disp = Show(gm_rim, view)
        gm_rim_disp.Representation = 'Wireframe'
        gm_rim_disp.ColorArrayName = ['POINTS', '']  # see ls_disp above for why
        gm_rim_color = GM_RIM_COLORS[panel_idx % len(GM_RIM_COLORS)]
        gm_rim_disp.AmbientColor = gm_rim_color
        gm_rim_disp.DiffuseColor = gm_rim_color
        gm_rim_disp.LineWidth = GM_RIM_LINE_WIDTH

        # Camera: transverse view, looking down -z (from further along the
        # scan direction back toward the cut face -- see header), up = -y so
        # atmosphere is "up" in the frame, matching the other two views'
        # convention. CameraParallelScale ties to the vertical (y) half-
        # extent, same rule the other two scripts use for their own vertical
        # screen axis.
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

    # Stack the panels side-by-side. All panels share one fixed crop window
    # (unlike the other two scripts), so every panel is already the same
    # size and fully populated -- no per-panel whitespace trim needed.
    import matplotlib.image as mpimg
    import numpy as np

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
