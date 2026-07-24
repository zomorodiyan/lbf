#!/usr/bin/env python3
"""Schematic of the shared VDEP domain, metal slab, laser track, and static
refinement boxes (2026-07-23 clipped geometry).

Values match report/jul23-tasks.md / tutorials/laserbeamFoam/vdep/testrun57_vdep_mesh_test.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({
    "font.size": 36,
    "axes.labelsize": 40,
    "legend.fontsize": 30,
    "xtick.labelsize": 30,
    "ytick.labelsize": 30,
})

# --- shared geometry (um), 2026-07-23 clipped domain ---
X0, X1 = -320, 320
Y0, Y1 = 0, 500            # was 0, 700
Z0, Z1 = 0, 3200
Y_SURF = 200                # was 300
TRACK_X = 0
TRACK_Z0, TRACK_Z1 = 500, 2900
LASER_R = 35

# static refinement boxes (see report/jul23-tasks.md / testrun57_vdep_mesh_test)
BOX1 = dict(x0=-100, x1=100, y0=200, y1=400, z0=700, z1=2700,
            color="forestgreen", label="refine box 1: 5µm, y[200,400]µm, x±100µm, z[700,2700]µm")
BOX2 = dict(x0=-250, x1=250, y0=145, y1=205, z0=Z0, z1=2700,
            color="royalblue", label="refine box 2 (surface overflow): 5µm, y[145,205]µm, x±250µm, z≤2700µm")

# laser positions to depict: (time label, z position)
Z_100US = TRACK_Z0 + 600  # 100 us * 6 um/us = 600 um
beam_events = [
    ("0 μs", TRACK_Z0),
    ("100 μs", Z_100US),
    ("400 μs", TRACK_Z1),
]

fig = plt.figure(figsize=(30, 22))
ax3d = fig.add_subplot(111, projection="3d")
# mplot3d's automatic depth-sort picks one centroid depth per artist and can draw a
# large translucent face after (on top of) a thin line that actually passes in front
# of it for most of its length -- the classic "line behind glass" artifact. Disabling
# it and setting explicit zorder per artist gives predictable, correct draw order for
# this kind of mostly-disjoint schematic (box shell behind, beams/labels in front).
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
                        linewidth=1.2, linestyle="-", label=None):
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

    first = True
    for i, j in ALL_EDGES:
        (ax0, ay0, az0), (ax1_, ay1_, az1_) = verts[i], verts[j]
        hidden = is_hidden(i, j)
        z = zorder_hidden if hidden else zorder_visible
        lbl = label if (first and not hidden) else None
        ax3d.plot([ax0, ax1_], [ay0, ay1_], [az0, az1_], color=color, linewidth=linewidth,
                   linestyle=linestyle, zorder=z, label=lbl)
        if lbl is not None:
            first = False


# domain wireframe -- visible edges on top (25), far-side edges behind the solid
# fills (1)
draw_box_wireframe(X0, X1, Y0, Y1, Z0, Z1, color="k", zorder_visible=25, zorder_hidden=1,
                    linewidth=1.6)

# metal slab (solid), y = Y_SURF -> Y1 (physical height). No edgecolors here: a
# Poly3DCollection draws all 6 of its own faces' edges without occlusion, which
# reproduced the same far-corner-lines-through-the-solid artifact as the domain
# box did. The domain wireframe above (which *is* occlusion-aware) already
# supplies the outer box lines; the metal/gas color change marks the surface.
metal = [[remap(*pt) for pt in face] for face in box_faces(X0, X1, Y_SURF, Y1, Z0, Z1)]
ax3d.add_collection3d(Poly3DCollection(metal, facecolors="silver", edgecolors="none",
                                        alpha=0.85, zorder=2))

# gas headspace, y = Y0 -> Y_SURF (physical height), drawn translucent for contrast
gas = [[remap(*pt) for pt in face] for face in box_faces(X0, X1, Y0, Y_SURF, Z0, Z1)]
ax3d.add_collection3d(Poly3DCollection(gas, facecolors="lightcyan", edgecolors="none",
                                        alpha=0.15, zorder=3))

# static refinement boxes -- drawn on top of the metal/gas fills (zorder > 3) but
# below the laser track/beams/labels, dashed so they read distinctly from the
# solid-black domain outline
draw_box_wireframe(**{k: v for k, v in BOX1.items() if k not in ("color", "label")},
                    color=BOX1["color"], zorder_visible=12, zorder_hidden=4,
                    linewidth=3.0, linestyle="--", label=BOX1["label"])
draw_box_wireframe(**{k: v for k, v in BOX2.items() if k not in ("color", "label")},
                    color=BOX2["color"], zorder_visible=13, zorder_hidden=4,
                    linewidth=3.0, linestyle="--", label=BOX2["label"])

# full laser track line along the surface (active portion, red)
tx, ty, tz = remap(TRACK_X, Y_SURF, TRACK_Z0)
tx2, ty2, tz2 = remap(TRACK_X, Y_SURF, TRACK_Z1)
ax3d.plot([tx, tx2], [ty, ty2], [tz, tz2], color="red", linewidth=3, zorder=15,
          label="laser track")

# track line extended beyond the active laser region to the domain edges
# (black, no legend entry) -- just shows where the surface/centerline continues
z0e, z0e2 = remap(TRACK_X, Y_SURF, Z0), remap(TRACK_X, Y_SURF, TRACK_Z0)
z1e, z1e2 = remap(TRACK_X, Y_SURF, TRACK_Z1), remap(TRACK_X, Y_SURF, Z1)
ax3d.plot([z0e[0], z0e2[0]], [z0e[1], z0e2[1]], [z0e[2], z0e2[2]],
          color="black", linewidth=2, zorder=14)
ax3d.plot([z1e[0], z1e2[0]], [z1e[1], z1e2[1]], [z1e[2], z1e2[2]],
          color="black", linewidth=2, zorder=14)

# laser beams (top of domain -> metal surface) at each timestamped location.
# A full plot_surface tube z-fights badly against the box's Poly3DCollections in
# mplot3d, so instead: a thick translucent centerline (renders like the other
# Line3D objects, e.g. the red track) plus a flat filled disk where it lands on
# the metal surface (coplanar with the metal face, so it sorts cleanly too).
theta = np.linspace(0, 2 * np.pi, 40)

for label, zc in beam_events:
    bx, by, bz = remap(TRACK_X, Y0, zc)
    sx, sy, sz = remap(TRACK_X, Y_SURF, zc)
    ax3d.plot([bx, sx], [by, sy], [bz, sz], color="orange", linewidth=6, alpha=0.85,
              solid_capstyle="round", zorder=16)

    disk = [remap(TRACK_X + LASER_R * np.cos(t), Y_SURF, zc + LASER_R * np.sin(t))
            for t in theta]
    ax3d.add_collection3d(Poly3DCollection([disk], facecolors="orangered",
                                            edgecolors="darkorange", linewidths=1.0,
                                            alpha=0.95, zorder=17))

    # time label, close above the top of the beam line
    lx, ly, lz = remap(TRACK_X, Y0 - 35, zc)
    ax3d.text(lx, ly, lz, label, color="darkorange", fontsize=34, fontweight="bold",
              ha="center", va="bottom", zorder=20)
    # (z-location is no longer labeled here -- the z-axis ticks now show 0.5/1.1/2.9
    # directly, which made this redundant)

ax3d.set_xlabel("x (µm)", labelpad=75)
ax3d.set_ylabel("z (µm)", labelpad=65)
ax3d.set_zlabel("y (µm)", labelpad=55)

# Ticks curated to only meaningful values (domain extent + feature locations),
# not an even grid -- e.g. x shows the domain half-width and both refinement
# boxes' half-widths; z shows the domain extent plus the three beam-event
# locations (which is also why the separate z-location text labels were
# dropped above).
X_TICKS = [X0, BOX2["x0"], BOX1["x0"], 0, BOX1["x1"], BOX2["x1"], X1]
Z_TICKS = sorted({Z0, BOX1["z0"], TRACK_Z0, Z_100US, BOX1["z1"], TRACK_Z1, Z1})
Y_TICKS = sorted({Y0, BOX2["y0"], BOX2["y1"], BOX1["y1"], Y1})  # BOX1["y0"] (200) omitted, redundant next to BOX2["y1"] (205)
def fmt(v):
    return f"{v:g}"


def stagger(values, extra_depth=None):
    # Alternate a leading blank line on every other label so dense tick runs
    # don't visually overlap, staggering them up/down instead of needing more
    # horizontal space. extra_depth: {index: n} adds n extra blank lines to
    # push specific labels further out -- used where the z=3200 and y=500
    # domain-edge labels meet at a shared corner and collide even after the
    # normal 1-line stagger.
    extra_depth = extra_depth or {}
    out = []
    for i, v in enumerate(values):
        depth = (1 if i % 2 else 0) + extra_depth.get(i, 0)
        out.append(("\n" * depth) + fmt(v))
    return out


ax3d.set_xticks(X_TICKS)
ax3d.set_xticklabels(stagger(X_TICKS))
ax3d.set_yticks(Z_TICKS)
ax3d.set_yticklabels([fmt(v) for v in Z_TICKS])
ax3d.set_zticks(Y_TICKS)
ax3d.set_zticklabels([fmt(v) for v in Y_TICKS])

# true mm ratios make this an unreadable sliver (scan extent is ~5x the other two) -- compress it
ax3d.set_box_aspect((X1 - X0, (Z1 - Z0) / 3.0, Y1 - Y0))
ax3d.invert_zaxis()  # y=0 (atmosphere) at top, y=Y1 (substrate) at bottom, visually
ax3d.xaxis.pane.set_facecolor("white")
ax3d.yaxis.pane.set_facecolor("white")
ax3d.zaxis.pane.set_facecolor("white")
ax3d.tick_params(axis="both", labelsize=38, pad=14)
ax3d.xaxis.set_tick_params(labelsize=34, pad=18)
ax3d.yaxis.set_tick_params(labelsize=34, pad=28)
ax3d.zaxis.set_tick_params(labelsize=38, pad=34)
ax3d.view_init(elev=20, azim=-55)
ax3d.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.02))

GRID_DEFAULT = "0.8"
CLOSE = 1e-9


def x_tick_color(v):
    if abs(v - BOX2["x0"]) < CLOSE or abs(v - BOX2["x1"]) < CLOSE:
        return BOX2["color"]
    if abs(v - BOX1["x0"]) < CLOSE or abs(v - BOX1["x1"]) < CLOSE:
        return BOX1["color"]
    return GRID_DEFAULT


def z_tick_color(v):
    if abs(v - BOX1["z0"]) < CLOSE or abs(v - BOX1["z1"]) < CLOSE:
        return BOX1["color"]
    if abs(v - TRACK_Z0) < CLOSE or abs(v - Z_100US) < CLOSE or abs(v - TRACK_Z1) < CLOSE:
        return "red"
    return GRID_DEFAULT


def y_tick_color(v):
    if abs(v - BOX2["y0"]) < CLOSE or abs(v - BOX2["y1"]) < CLOSE:
        return BOX2["color"]
    if abs(v - BOX1["y0"]) < CLOSE or abs(v - BOX1["y1"]) < CLOSE:
        return BOX1["color"]
    return GRID_DEFAULT


def nudge_label(lab, dx_chars, dy_chars, fontsize):
    """Permanently shift a 3D tick label by a number of character-widths in
    screen space (right/up positive). A plain set_transform() offset gets
    silently discarded: mplot3d's Axis3D.draw() resets each tick label's
    transform back to plain transData on every draw (including the one
    savefig() triggers internally). Patching set_transform so it composes
    our offset on top of *whatever* transform it's given next survives that,
    since it re-applies on every future call, not just once now."""
    char_px = 0.6 * fontsize
    nudge = ScaledTranslation(dx_chars * char_px / 72, dy_chars * char_px / 72,
                               fig.dpi_scale_trans)
    orig_set_transform = lab.set_transform

    def _patched(t, orig=orig_set_transform, nudge=nudge):
        orig(t + nudge)

    lab.set_transform = _patched
    lab.set_transform(lab.get_transform())  # apply immediately too


# Color the box-boundary x/y/z tick labels to match their box's outline color
fig.canvas.draw()
for loc, lab in zip(ax3d.xaxis.get_majorticklocs(), ax3d.xaxis.get_ticklabels()):
    c = x_tick_color(loc)
    if c != GRID_DEFAULT:
        lab.set_color(c)
    if abs(loc - BOX2["x0"]) < CLOSE:
        nudge_label(lab, dx_chars=-2, dy_chars=0, fontsize=34)
    elif abs(loc - BOX2["x1"]) < CLOSE:
        nudge_label(lab, dx_chars=-2, dy_chars=1, fontsize=34)
for loc, lab in zip(ax3d.yaxis.get_majorticklocs(), ax3d.yaxis.get_ticklabels()):
    c = z_tick_color(loc)
    if c != GRID_DEFAULT:
        lab.set_color(c)
    nudge_label(lab, dx_chars=1, dy_chars=1, fontsize=34)
    if abs(loc - BOX1["z0"]) < CLOSE:  # 700: extra nudge, up+right
        nudge_label(lab, dx_chars=0.5, dy_chars=0.5, fontsize=34)
    elif abs(loc - BOX1["z1"]) < CLOSE:  # 2700: extra nudge, down+left
        nudge_label(lab, dx_chars=-0.5, dy_chars=-0.5, fontsize=34)
for loc, lab in zip(ax3d.zaxis.get_majorticklocs(), ax3d.zaxis.get_ticklabels()):
    c = y_tick_color(loc)
    if c != GRID_DEFAULT:
        lab.set_color(c)

fig.tight_layout()
fig.savefig("report/domain_schematic.png", dpi=160)
print("wrote report/domain_schematic.png")
