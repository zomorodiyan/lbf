#!/usr/bin/env python3
"""Schematic of the shared VDEP domain, metal slab, and laser track/position.

Values match report/jul23-tasks.md / tutorials/laserbeamFoam/vdep/testrun57_vdep_mesh_test.

Usage: plot_domain_schematic.py [laser_z_um] [time_label] [output_png]
All three are optional (defaults below) so this still runs standalone for a
quick look. Called per-frame by render_view.py's render_transverse(), which
passes the *actual* current laser z position and a time label for that
frame -- the single beam/cylinder marker always reflects the specific
timestep being rendered, not a fixed set of illustrative timestamps.
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
TIME_LABEL = sys.argv[2] if len(sys.argv) > 2 else "0 us"
OUTPUT_PNG = sys.argv[3] if len(sys.argv) > 3 else "results/domain_schematic.png"

# The one deliberately-large font size (user request, 2026-08-02: "font's
# much larger"), reused for the x/y/z axis labels too ("use the same for
# the other values in the plot") since axis tick *numbers* are dropped
# entirely (see set_xticks([]) etc. below), leaving the time label and axis
# labels as the only text in the schematic.
TIME_FONTSIZE = 70

fig = plt.figure(figsize=(30, 22))
ax3d = fig.add_subplot(111, projection="3d")
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
ax3d.add_collection3d(Poly3DCollection(metal, facecolors="gray", edgecolors="none",
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

# Cross-section cut-line markers, matching render_view.py's transverse
# panels (OFFSETS_BEHIND_LASER behind the laser, colored via
# GM_RIM_COLORS) -- duplicated here rather than imported (same "small
# constants duplicated across each standalone pvpython script" convention
# already noted in render_view.py's own top_screenshot comment; keep in
# sync by hand if either changes). Each cut is drawn as an L-shape: along
# the plate's top surface (spanning the full x width) and down the visible
# x=X1 side face (spanning the plate's own solid depth, Y_SURF->Y1) at
# that cut's z -- user request, 2026-08-02 ("add cross-section lines to
# the schematic ... on the plate's top surface and the x=0.32 surface of
# the plate"). x=X1 (not X0) because that's this view's visible/front
# face -- draw_box_wireframe's own is_hidden() already treats x=X0 as the
# hidden far side for this view angle.
OFFSETS_BEHIND_LASER_UM = [700, 600, 500]
GM_RIM_COLORS = [(0.0, 1.0, 1.0), (0.0, 0.8, 0.0), (1.0, 1.0, 0.0)]
for _offset_um, _color in zip(OFFSETS_BEHIND_LASER_UM, GM_RIM_COLORS):
    _z_cut = LASER_Z_UM - _offset_um
    if not (Z0 <= _z_cut <= Z1):
        continue
    _top0 = remap(X0, Y_SURF, _z_cut)
    _top1 = remap(X1, Y_SURF, _z_cut)
    ax3d.plot([_top0[0], _top1[0]], [_top0[1], _top1[1]], [_top0[2], _top1[2]],
              color=_color, linewidth=4, zorder=15)
    _side0 = remap(X1, Y_SURF, _z_cut)
    _side1 = remap(X1, Y1, _z_cut)
    ax3d.plot([_side0[0], _side1[0]], [_side0[1], _side1[1]], [_side0[2], _side1[2]],
              color=_color, linewidth=4, zorder=15)

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
BOX_ASPECT = (X1 - X0, (Z1 - Z0) / 3.0, Y1 - Y0)
if hasattr(ax3d, 'set_box_aspect'):  # matplotlib >=3.3 only -- this Docker
    ax3d.set_box_aspect(BOX_ASPECT)  # image's pvpython has 3.1.1, doesn't have this
else:
    # Pre-3.3 workaround: without set_box_aspect, mplot3d maps each axis's
    # own view-limit *span* to the same on-screen cube edge regardless of
    # its actual data range -- so all 3 axes default to looking equal-sized
    # (confirmed by testing: the box rendered as a plain cube here, visibly
    # different from the intended elongated-in-z proportions -- user-
    # flagged, 2026-08-02, "don't change the dimensions of the cube").
    # Faking the same box_aspect ratios is possible by *padding* each
    # axis's view limits (not touching any plotted data) beyond its true
    # data range: an axis padded to span S instead of its true range T
    # occupies only T/S of its own cube edge, so choosing per-axis padding
    # ratios that match BOX_ASPECT reproduces the same relative on-screen
    # proportions. The axis with the single largest BOX_ASPECT value needs
    # zero padding (it already IS the longest edge); the other two get
    # padded up to match it. Derivation: for target ratio R_i and true
    # range T_i, span S_i = T_i * max(R)/R_i keeps (T_i/S_i) proportional
    # to R_i for every axis, and is always >= T_i (a valid, non-clipping
    # span) since max(R)/R_i >= 1.
    true_ranges = (X1 - X0, Z1 - Z0, Y1 - Y0)
    centers = ((X0 + X1) / 2, (Z0 + Z1) / 2, (Y0 + Y1) / 2)
    max_r = max(BOX_ASPECT)
    spans = [t * max_r / r for t, r in zip(true_ranges, BOX_ASPECT)]
    (sx0, sx1), (sy0, sy1), (sz0, sz1) = (
        (c - s / 2, c + s / 2) for c, s in zip(centers, spans))
    ax3d.set_xlim3d(sx0, sx1)
    ax3d.set_ylim3d(sy0, sy1)
    ax3d.set_zlim3d(sz0, sz1)
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

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=160)
print(f"wrote {OUTPUT_PNG}")
