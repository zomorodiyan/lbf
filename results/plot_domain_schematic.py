#!/usr/bin/env python3
"""Schematic of the shared VDEP domain, metal slab, and laser track/position.

Values match report/jul23-tasks.md / tutorials/laserbeamFoam/vdep/testrun57_vdep_mesh_test.

Usage: plot_domain_schematic.py [laser_z_um] [time_us] [output_png]
All three are optional (defaults below) so this still runs standalone for a
quick look. Called per-frame by render_view.py's render_transverse(), which
passes the *actual* current laser z position and elapsed time (in
microseconds, as a plain number) for that frame -- the single beam/cylinder
marker always reflects the specific timestep being rendered, not a fixed
set of illustrative timestamps.

time_us is a bare number, not a pre-formatted "N μs" string -- this script
appends its own "μs" suffix (below) rather than receiving one through argv,
since this Docker image's pvpython (Python 3.6) encodes subprocess argv
using the *caller's* own locale (ASCII/POSIX here), which raises
UnicodeEncodeError on a μ character before the child process ever starts;
building the suffix here, as a literal in this script's own source, sidesteps
that entirely (user report, 2026-08-03).
"""
import sys

import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({
    "font.size": 36,
})

# --- shared geometry (um), 2026-07-23 clipped domain ---
X0, X1 = -320, 320
Y0, Y1 = 0, 500            # was 0, 700
Z0, Z1 = 0, 3200
Y_SURF = 200                # was 300
TRACK_X = 0
TRACK_Z0, TRACK_Z1 = 500, 2900

# Laser position to depict -- one cylinder only, at the actual timestep
# being rendered (see module docstring). CLI overrides; defaults below let
# this still run standalone.
LASER_Z_UM = float(sys.argv[1]) if len(sys.argv) > 1 else TRACK_Z0
TIME_US = sys.argv[2] if len(sys.argv) > 2 else "0"
TIME_LABEL = f"{TIME_US} μs"  # was "us" -- user request, 2026-08-03
OUTPUT_PNG = sys.argv[3] if len(sys.argv) > 3 else "results/domain_schematic.png"

# The one deliberately-large font size (user request, 2026-08-02: "font's
# much larger"), reused for the x/y/z axis labels too ("use the same for
# the other values in the plot") since axis tick *numbers* are dropped
# entirely (see set_xticks([]) etc. below), leaving the time label and axis
# labels as the only text in the schematic.
TIME_FONTSIZE = 70

fig = plt.figure(figsize=(30, 22))
# add_axes([0,0,1,1]), not add_subplot(111) -- occupies the *entire* figure
# canvas with no default matplotlib subplot margins, so this axes' own
# pixel-box aspect ratio is exactly the figure's (30:22), a fixed, known
# quantity. That matters because mplot3d stretches its rendered content
# independently to fill whatever pixel box the axes occupies (no equal-
# aspect enforcement), so a given 3D direction (e.g. the z/scan axis) only
# renders at the *same on-screen slope* in two separate figures if their
# axes boxes share the same pixel aspect ratio -- true here and in
# proto_transverse_3d.py's own fig.add_axes([0,0,1,1]) (matched figsize
# ratio, same reasoning), by user request, 2026-08-03: "modify the
# schematic camera position ... so it matches that of the cross-section
# plot ... I want them to be parallel". add_subplot(111)'s default margins
# (further perturbed here by the tight_layout() call this replaces) were
# never a reliably matching, known ratio.
ax3d = fig.add_axes([0, 0, 1, 1], projection="3d")
# mplot3d's automatic depth-sort picks one centroid depth per artist and can draw a
# large translucent face after (on top of) a thin line that actually passes in front
# of it for most of its length -- the classic "line behind glass" artifact. Disabling
# it and setting explicit zorder per artist gives predictable, correct draw order for
# this kind of mostly-disjoint schematic (box shell behind, beam/label in front).
ax3d.computed_zorder = False


def remap(x, y, z):
    # mplot3d always draws its 3rd argument as the vertical screen axis; physical y
    # (height/depth) should read as "up", so remap (x, y, z) -> (x, z, y) for plotting.
    return x, z, y


def box_faces(x0, x1, y0, y1, z0, z1):
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    return [
        [v[0], v[1], v[2], v[3]],  # z0
        [v[4], v[5], v[6], v[7]],  # z1
        [v[0], v[1], v[5], v[4]],  # y0
        [v[3], v[2], v[6], v[7]],  # y1
        [v[0], v[3], v[7], v[4]],  # x0
        [v[1], v[2], v[6], v[5]],  # x1
    ]


def box_vertices(x0, x1, y0, y1, z0, z1):
    return [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]


ALL_EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]


def apply_box_aspect(ax, x0, x1, y0, y1, z0, z1, aspect):
    """Set `ax`'s (already-plotted, i.e. x/z/y-remapped) view limits so its
    3 axes render in the given relative proportions -- shared by the main
    domain axes and the zoomed-inset axes below (both need this, just with
    different true_ranges/aspect). See set_box_aspect's own call site
    (below) for why the pre-3.3 workaround branch exists at all."""
    if hasattr(ax, 'set_box_aspect'):  # matplotlib >=3.3 only -- this
        ax.set_xlim3d(x0, x1)          # Docker image's pvpython has 3.1.1,
        ax.set_ylim3d(y0, y1)          # doesn't have this
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


def draw_box_wireframe(x0, x1, y0, y1, z0, z1, color, zorder_visible, zorder_hidden,
                        linewidth=1.2, linestyle="-"):
    """Draw a box's 12 edges, sending the far-side ones (on its own x=x0 or z=z1
    face) to a lower zorder so a solid object in front of them still occludes
    correctly -- same rule established for the outer domain box, for this fixed
    view (elev=20, azim=-55)."""
    raw_verts = box_vertices(x0, x1, y0, y1, z0, z1)
    verts = [remap(*pt) for pt in raw_verts]

    def is_hidden(i, j):
        (px0, _, pz0) = raw_verts[i]
        (px1, _, pz1) = raw_verts[j]
        return (px0 == x0 and px1 == x0) or (pz0 == z1 and pz1 == z1)

    for i, j in ALL_EDGES:
        (ax0, ay0, az0), (ax1_, ay1_, az1_) = verts[i], verts[j]
        hidden = is_hidden(i, j)
        z = zorder_hidden if hidden else zorder_visible
        ax3d.plot([ax0, ax1_], [ay0, ay1_], [az0, az1_], color=color, linewidth=linewidth,
                   linestyle=linestyle, zorder=z)


# domain wireframe -- visible edges on top (25), far-side edges behind the solid
# fills (1)
draw_box_wireframe(X0, X1, Y0, Y1, Z0, Z1, color="k", zorder_visible=25, zorder_hidden=1,
                    linewidth=5.0)  # thickened 1.6 -> 2.5 -> 5.0 (2x further,
                                    # user request, 2026-08-02)

# metal slab (solid), y = Y_SURF -> Y1 (physical height). No edgecolors here: a
# Poly3DCollection draws all 6 of its own faces' edges without occlusion, which
# reproduced the same far-corner-lines-through-the-solid artifact as the domain
# box did. The domain wireframe above (which *is* occlusion-aware) already
# supplies the outer box lines (the plate's own bottom/side edges, shared
# with the domain box at Y1 and at X0/X1/Z0/Z1) -- but not the plate's own
# *top* rim at y=Y_SURF, partway up the domain box, which used to be
# implied only by the color contrast against the (now-removed) translucent
# gas-headspace fill above it. Without that fill the rim had no line at
# all -- "the edge of the plate is non-existent" (user report, 2026-08-02)
# -- so it's drawn explicitly below.
metal = [[remap(*pt) for pt in face] for face in box_faces(X0, X1, Y_SURF, Y1, Z0, Z1)]
# Light gray, was plain "gray" -- darker gray made the cross-section
# outline rectangles (drawn on top, see the inset below) hard to pick out
# against it (user request, 2026-08-03).
ax3d.add_collection3d(Poly3DCollection(metal, facecolors="0.88", edgecolors="none",
                                        alpha=0.85, zorder=2))

# Plate top-rim outline at y=Y_SURF -- draw_box_wireframe with y0=y1=Y_SURF
# collapses box_vertices() into a degenerate (zero-height) box: the 4
# "vertical" edges connecting its (identical) top/bottom vertex pairs
# become zero-length and invisible, leaving exactly the 4 perimeter edges
# of the y=Y_SURF rectangle -- each still run through is_hidden()'s own
# occlusion logic, so the near/right edges render visible and the
# far/left ones render behind the solid, same as every other edge in this
# figure. Cheaper than writing a separate one-off 4-edge helper. Dark gray
# (not black, unlike the outer domain wireframe) -- this rim is a plate
# edge, not an actual domain-boundary edge, and should read as visually
# distinct from the ones that are (user request, 2026-08-02; "dimgray"
# darkened further to (0.25,0.25,0.25) per user follow-up).
draw_box_wireframe(X0, X1, Y_SURF, Y_SURF, Z0, Z1, color=(0.25, 0.25, 0.25), zorder_visible=25, zorder_hidden=1,
                    linewidth=5.0)

# Cross-section window markers, matching render_view.py's transverse overlay
# (X_CROSS_SECTIONS, colored via CROSS_SECTION_COLORS) -- duplicated here
# rather than imported (same "small constants duplicated across each
# standalone pvpython script" convention already noted in render_view.py's
# own top_screenshot comment; keep in sync by hand if either changes).
# Redesigned 2026-08-03 (from the old distance-behind-laser z-cut
# "L-shape" markers, spanning the full plate) to match render_transverse's
# X-based cross-sections: each is now a small rectangular window (not
# a full-width/full-depth line) at a fixed x, spanning the shared y/z crop
# window used by every cross-section -- drawn as a wireframe rectangle
# lying in the plane normal to x at that x position. Reduced from 5 back
# to 3 positions same-day (user request, 2026-08-03) -- the middle 3 of
# the earlier 5-position set (+-250um came back empty at every timestep
# tested: the surface never reaches that far off-center within the shared
# y=[250,350]um window). Z is a *fixed* absolute window
# (Z_WINDOW_MIN_UM/MAX_UM), not laser-relative, so these markers stay put
# across an entire batch of frames regardless of the laser's own (still
# separately drawn) position.
X_CROSS_SECTIONS_UM = [-25, 0, 25]
Z_WINDOW_MIN_UM, Z_WINDOW_MAX_UM = 1500, 2000
Y_MIN_UM, Y_MAX_UM = 225, 350
CROSS_SECTION_COLORS = [
    (0.85, 0.10, 0.10),  # -25um red
    (0.10, 0.60, 0.10),  #   0um green
    (0.10, 0.45, 0.90),  # +25um blue
]
MARKER_PTS_3D = []  # every marker corner (remapped) -- the tiny markers
                     # used to be drawn directly in the main scene here, but
                     # at the domain's own scale they were an illegible
                     # smudge (user report, 2026-08-03, see Media.jpg: a
                     # hand-drawn sketch of a callout box zooming in on this
                     # cluster instead). Per user follow-up, the tiny
                     # in-place markers are no longer drawn at all -- only
                     # their corner points are kept, to place the callout
                     # rectangle and zoomed inset below at the right spot.
for _x_cut, _color in zip(X_CROSS_SECTIONS_UM, CROSS_SECTION_COLORS):
    _z0 = max(Z0, min(Z1, Z_WINDOW_MIN_UM))
    _z1 = max(Z0, min(Z1, Z_WINDOW_MAX_UM))
    if _z0 >= _z1:
        continue  # window not within the domain (shouldn't happen; fixed values)
    _corners = [
        (_x_cut, Y_MIN_UM, _z0), (_x_cut, Y_MIN_UM, _z1),
        (_x_cut, Y_MAX_UM, _z1), (_x_cut, Y_MAX_UM, _z0),
        (_x_cut, Y_MIN_UM, _z0),
    ]
    MARKER_PTS_3D.extend(remap(*pt) for pt in _corners)

# Laser track line along the surface (active portion) -- no longer
# extended to the domain edges in black (removed, user request,
# 2026-08-02: "no need to extend the laser track all the way"). Split at
# the laser's current position: red for the already-traveled portion
# (TRACK_Z0 to the laser), dotted white for the remaining portion ahead of
# it (up to TRACK_Z1) -- user request, 2026-08-02 ("red up to the laser
# location and dotted white line for the rest of the track"). Clamped to
# [TRACK_Z0,TRACK_Z1] in case the laser hasn't started yet or has already
# passed the end of the active region.
_laser_track_z = max(TRACK_Z0, min(TRACK_Z1, LASER_Z_UM))
tx, ty, tz = remap(TRACK_X, Y_SURF, TRACK_Z0)
tx2, ty2, tz2 = remap(TRACK_X, Y_SURF, _laser_track_z)
ax3d.plot([tx, tx2], [ty, ty2], [tz, tz2], color="red", linewidth=4.5, zorder=15)  # thickened (was 3)
if _laser_track_z < TRACK_Z1:
    tx3, ty3, tz3 = remap(TRACK_X, Y_SURF, TRACK_Z1)
    ax3d.plot([tx2, tx3], [ty2, ty3], [tz2, tz3], color="white", linewidth=4.5,
              linestyle=':', zorder=15)

# Single laser beam (top of domain -> metal surface) at the actual current
# position (LASER_Z_UM): a thick translucent centerline (renders like the
# other Line3D objects, e.g. the red track) plus a marker dot where it lands
# on the metal surface. The landing marker used to be a small flat disk at
# the beam's true physical radius (LASER_R=35um) -- but without
# set_box_aspect (matplotlib 3.1.1 in this Docker image doesn't have it,
# see below) the z-axis (scan direction, 3200um) is squeezed down to fit
# the same visual cube edge as x (640um) instead of the previously-tuned
# ~1.4x-wider proportion, shrinking that disk to near-invisible. A
# fixed-pixel-size marker sidesteps that entirely, since its size doesn't
# depend on the surrounding (now more compressed) data scale.
bx, by, bz = remap(TRACK_X, Y0, LASER_Z_UM)
sx, sy, sz = remap(TRACK_X, Y_SURF, LASER_Z_UM)
ax3d.plot([bx, sx], [by, sy], [bz, sz], color="orange", linewidth=10, alpha=1.0,
          solid_capstyle="round", zorder=16)  # fully solid + thicker (user
                                               # request, 2026-08-02) -- was
                                               # alpha=0.85, linewidth=6
ax3d.plot([sx], [sy], [sz], marker="o", markersize=30, markerfacecolor="orangered",
          markeredgecolor="darkorange", markeredgewidth=2, zorder=17)

# Time label, above the top of the beam line -- much larger than the
# rest of the schematic's text (user request, 2026-08-02: "font's much
# larger"). Offset raised 35 -> 70 -> 100um, progressively further from
# the beam's own top (user follow-up requests, 2026-08-02: "move the time
# number a little bit further from laser cylinder", then "a bit higher").
# Color switched orange -> black and size bumped past TIME_FONTSIZE
# (user request, 2026-08-02: "show it in black color and slightly larger"),
# then multiplied by another 1.5x on top of that (user request,
# 2026-08-02: "increase the time font size to 1.5 times of the current
# value") -- net multiplier on TIME_FONTSIZE is 1.15 * 1.5 = 1.725.
lx, ly, lz = remap(TRACK_X, Y0 - 100, LASER_Z_UM)
ax3d.text(lx, ly, lz, TIME_LABEL, color="black", fontsize=round(TIME_FONTSIZE * 1.725),
          fontweight="bold", ha="center", va="bottom", zorder=20)

# x/y/z labels placed by hand, near the midpoint of one representative edge
# per axis and pushed a bit further outward from the box's own center --
# user request, 2026-08-02 ("very close to the domain boundaries at
# appropriate location"). set_xlabel/set_ylabel/set_zlabel (with generous
# labelpad) left them floating noticeably far from the box and didn't
# track the rotation predictably (same class of view-angle-dependent
# mplot3d placement logic noted in the dropped arrow attempt above); manual
# placement guarantees they land where intended regardless of view angle.
_box_center = ((X0 + X1) / 2, (Y0 + Y1) / 2, (Z0 + Z1) / 2)
LABEL_PUSH_UM = 60  # extra outward push past the edge midpoint
# Bigger than TIME_FONTSIZE, and each label carries a direction arrow --
# user request, 2026-08-02 ("increase the font size of x, y, z", "change
# them to x |> y ^ z /^", then a follow-up "increase the font sizes ...
# further" bumped the multiplier 1.3 -> 1.6). Arrow glyphs are plain
# Unicode (DejaVu Sans, this Docker image's default font, covers them --
# unlike the missing micro-sign glyph noted elsewhere in this repo).
AXIS_FONTSIZE = round(TIME_FONTSIZE * 1.6)
AXIS_LABELS = {"x": "x↘", "y": "y↓", "z": "z↗"}
for axis_name, edge_mid in [
    ("x", ((X0 + X1) / 2, Y1, Z0)),
    ("y", (X0, (Y0 + Y1) / 2, Z0)),
    ("z", (X0, Y1, (Z0 + Z1) / 2)),
]:
    dx, dy, dz = (edge_mid[i] - _box_center[i] for i in range(3))
    norm = (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5
    pushed = tuple(edge_mid[i] + LABEL_PUSH_UM * (edge_mid[i] - _box_center[i]) / norm for i in range(3))
    lx, ly, lz = remap(*pushed)
    label = ax3d.text(lx, ly, lz, AXIS_LABELS[axis_name], color="black", fontsize=AXIS_FONTSIZE,
                       fontweight="bold", ha="center", va="center", zorder=27)
    # Per-axis screen-space nudges (character units, not 3D data
    # coordinates) -- deliberately screen-space so the offset stays fixed
    # on the page regardless of any later view_init/azim change ("if there
    # was not rotation" -- user request, 2026-08-02), same ScaledTranslation
    # technique the (now-removed) tick-label nudging used, minus the
    # "reapply every draw" patch that was only needed for mplot3d's own
    # managed tick Text objects (reset by Axis3D.draw() every frame); this
    # is a plain ax3d.text() artist, unaffected by that.
    char_px = 0.6 * TIME_FONTSIZE
    dx_chars, dy_chars = {
        "x": (0, -1),   # 1 character down
        "y": (-1, 0),   # 1 character to the left
        "z": (12, -4),  # 12 characters to the right (5, then +4, then +3
                        # more per user follow-ups), plus the existing
                        # 4-characters-down nudge from the previous round
    }[axis_name]
    offset = ScaledTranslation(dx_chars * char_px / 72, dy_chars * char_px / 72, fig.dpi_scale_trans)
    label.set_transform(label.get_transform() + offset)

# No tick numbers on any axis -- user request, 2026-08-02 ("no need for
# axis to have values on them"). The axis labels above are enough context;
# this also means the box-aspect workaround below no longer needs to relabel
# a compressed z-axis with its true values.
ax3d.set_xticks([])
ax3d.set_yticks([])
ax3d.set_zticks([])
# mplot3d still draws its own axis spine line for each axis even with zero
# ticks (independent of our hand-drawn domain wireframe box) -- redundant
# now that nothing else aligns to it, and it can visibly jut out past our
# box since it's positioned using the padded set_xlim3d/etc. limits above,
# not the true domain extent. Hide it -- just the labels are wanted (user
# request, 2026-08-02: "no need to include x, y, z lines, just the labels").
ax3d.xaxis.line.set_visible(False)
ax3d.yaxis.line.set_visible(False)
ax3d.zaxis.line.set_visible(False)

# true mm ratios make this an unreadable sliver (scan extent is ~5x the other two) -- compress it
# (Pre-3.3 fallback: without set_box_aspect, mplot3d maps each axis's own
# view-limit *span* to the same on-screen cube edge regardless of its
# actual data range -- so all 3 axes default to looking equal-sized
# unless the view limits are padded to fake the target ratios -- see
# apply_box_aspect()'s own docstring.)
apply_box_aspect(ax3d, X0, X1, Z0, Z1, Y0, Y1, (X1 - X0, (Z1 - Z0) / 3.0, Y1 - Y0))
ax3d.invert_zaxis()  # y=0 (atmosphere) at top, y=Y1 (substrate) at bottom, visually
ax3d.xaxis.pane.set_facecolor("white")
ax3d.yaxis.pane.set_facecolor("white")
ax3d.zaxis.pane.set_facecolor("white")
ax3d.view_init(elev=20, azim=-40)  # rotated around the (physical) y axis:
                                    # -55 (original) -> -20 ("more left to
                                    # right") -> -30 -> -35 -> -40 (another
                                    # slight rotation, same direction as the
                                    # last two steps, user request,
                                    # 2026-08-02)

# Zoomed inset for the cross-section markers -- user request, 2026-08-03,
# after sharing a hand-drawn sketch (Media.jpg) of a callout box zooming in
# on this cluster, since at the domain's own scale the 3 markers are only
# 50um apart across a 640um-wide, 3200um-long box and read as a tiny,
# illegible smudge. A separate small 3D axes placed over that footprint
# (found by projecting the markers' 3D corners through the now-final camera
# transform), showing just the 3 cross-section planes as thin colored
# outlines -- no fill/hatch and no border rectangle around the inset itself
# (both dropped, user request, 2026-08-03: "simply colored outlines",
# "remove the purple square") -- at a scale large enough to make out their
# relative offset.
from mpl_toolkits.mplot3d import proj3d  # noqa: E402 (local import -- only
                                          # needed here, at the bottom, once
                                          # the main axes' projection is final)

fig.canvas.draw()  # force a draw so transData/get_proj reflect the final
                    # axes position (fixed at [0,0,1,1] now, but the camera
                    # projection itself still needs a draw to settle)


def project_to_fig_frac(pt3d):
    """3D data coords (already x/z/y-remapped) -> figure-fraction (0-1)
    coords, via the axes' current camera projection then its 2D transData."""
    x2, y2, _ = proj3d.proj_transform(*pt3d, ax3d.get_proj())
    disp = ax3d.transData.transform((x2, y2))
    return fig.transFigure.inverted().transform(disp)


_fig_pts = [project_to_fig_frac(pt) for pt in MARKER_PTS_3D]
_fxs, _fys = zip(*_fig_pts)
_marker_cx, _marker_cy = (min(_fxs) + max(_fxs)) / 2, (min(_fys) + max(_fys)) / 2

# The callout/inset box is centered on the markers' own on-screen position
# but drawn much bigger (a "magnifying glass held over the spot", per the
# sketch) rather than tucked in an empty corner -- CALLOUT_W chosen to
# comfortably fit 3 labeled swatch-sized planes while still clearly reading
# as "zoomed in on that spot" rather than an unrelated separate panel.
#
# CALLOUT_H is *derived*, not a second free constant: a 3D Axes stretches
# its own normalized projection independently to fill whatever *pixel* box
# it occupies (no equal-aspect enforcement here), so two Axes3D objects
# with the same view_init but *different* pixel-box aspect ratios render
# the same 3D direction at two different on-screen slopes -- per-axis data
# scaling (box_aspect/apply_box_aspect) never affects this, only the pixel
# box shape does. That was why the inset's cross-section rectangles didn't
# look parallel to the main scene's z-direction laser track (user report,
# 2026-08-03): CALLOUT_W/CALLOUT_H (0.30/0.36) didn't match ax3d's own
# actual (post-tight_layout) box aspect. Matching it here -- computed from
# ax3d.get_position() in inches, not hardcoded -- makes the two Axes3D
# objects apply the *same* pixel stretch, so a given 3D direction (e.g.
# physical z) lands at the same screen slope in both.
_main_pos = ax3d.get_position()
_fig_w_in, _fig_h_in = fig.get_size_inches()
_main_aspect = (_main_pos.width * _fig_w_in) / (_main_pos.height * _fig_h_in)
CALLOUT_W = 0.30
CALLOUT_H = (CALLOUT_W * _fig_w_in) / (_main_aspect * _fig_h_in)
callout_x0 = min(max(_marker_cx - CALLOUT_W / 2, 0.0), 1.0 - CALLOUT_W)
callout_y0 = min(max(_marker_cy - CALLOUT_H / 2, 0.0), 1.0 - CALLOUT_H)

ax_inset = fig.add_axes([callout_x0, callout_y0, CALLOUT_W, CALLOUT_H], projection="3d")
ax_inset.computed_zorder = False
ax_inset.patch.set_alpha(0)  # transparent background -- shows the main
                              # scene through any part of the inset's own
                              # (rectangular) bounding box not covered by
                              # its 3D content, so it reads as an overlay,
                              # not an opaque separate panel

INSET_PAD_X, INSET_PAD_YZ = 15, 25  # um, headroom around the 3 planes'
                                     # own bounding box so they don't touch
                                     # the inset's own edge
_ix0, _ix1 = min(X_CROSS_SECTIONS_UM) - INSET_PAD_X, max(X_CROSS_SECTIONS_UM) + INSET_PAD_X
_iy0, _iy1 = Y_MIN_UM - INSET_PAD_YZ, Y_MAX_UM + INSET_PAD_YZ
_iz0, _iz1 = Z_WINDOW_MIN_UM - INSET_PAD_YZ, Z_WINDOW_MAX_UM + INSET_PAD_YZ

# Only the middle/green (x=0) cross-section frame -- user request,
# 2026-08-03 -- at the same thickness as the laser track line (linewidth=
# 4.5, see its own ax3d.plot call above) rather than this inset's old
# linewidth=3.
_x_cut, _color = 0.0, CROSS_SECTION_COLORS[1]
_corners = [
    remap(_x_cut, Y_MIN_UM, Z_WINDOW_MIN_UM),
    remap(_x_cut, Y_MIN_UM, Z_WINDOW_MAX_UM),
    remap(_x_cut, Y_MAX_UM, Z_WINDOW_MAX_UM),
    remap(_x_cut, Y_MAX_UM, Z_WINDOW_MIN_UM),
    remap(_x_cut, Y_MIN_UM, Z_WINDOW_MIN_UM),
]
_xs, _ys, _zs = zip(*_corners)
ax_inset.plot(_xs, _ys, _zs, color=_color, linewidth=4.5)

# True (uncompressed) local aspect ratios -- unlike the main box, this
# region isn't dominated by the long scan axis, so no /3.0 fudge is needed;
# the 3 planes' actual relative proportions are what we want visible here.
# (Per-axis data scaling never changes an axis-aligned vector's *direction*,
# only its length, so this ratio has no bearing on whether the inset's
# z-direction edges look parallel to the main scene's laser track -- see
# CALLOUT_W/CALLOUT_H below for what actually controls that.)
apply_box_aspect(ax_inset, _ix0, _ix1, _iz0, _iz1, _iy0, _iy1,
                  (_ix1 - _ix0, _iz1 - _iz0, _iy1 - _iy0))
ax_inset.invert_zaxis()
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax_inset.set_zticks([])
ax_inset.xaxis.line.set_visible(False)
ax_inset.yaxis.line.set_visible(False)
ax_inset.zaxis.line.set_visible(False)
ax_inset.xaxis.pane.set_visible(False)
ax_inset.yaxis.pane.set_visible(False)
ax_inset.zaxis.pane.set_visible(False)
ax_inset.view_init(elev=20, azim=-40)  # same tilt as the main scene, so the
                                        # zoomed planes read as the same
                                        # 3D objects, just enlarged

fig.savefig(OUTPUT_PNG, dpi=160)
print(f"wrote {OUTPUT_PNG}")
