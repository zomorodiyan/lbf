# Lateral (through-thickness) X-ray-style projection of the melt pool --
# liquid/solid variant.
#
# Same surface-based ray-tracing approach as lateral_xray.py (see that file
# for the full history: sawtooth-elimination via marching-cubes isosurfaces
# + exact ray/surface intersection, and the vtkStaticCellLocator locator
# gotchas). This version additionally traces the liquid/solid (epsilon1)
# boundary and draws it as a dotted line over the same grayscale image --
# "the bottom of the melt pool" as seen from this lateral view -- without
# baking the liquid/solid distinction into the displayed grayscale itself
# (MU_SOLID == MU_LIQUID for the main image; see below for why). The
# gas/metal (vapor depression) boundary line that lateral_xray.py draws was
# dropped here by request -- only the melt-pool line remains.
#
# History of the boundary-line method, since this went through a few wrong
# turns:
#   v1: baked the liquid/solid difference directly into the *displayed*
#       attenuation (distinct MU_SOLID/MU_LIQUID), so the melt pool showed up
#       as a different shade with no separate line. Reverted: physically,
#       solid and liquid Al6061 are nearly the same density (the solver
#       treats them as one incompressible phase, see transportProperties'
#       single rho=2700 for "metal" -- only a ~0.5% thermal-expansion
#       correction distinguishes them), and baking the distinction into the
#       displayed image let ls_poly's degenerate axis-aligned slivers (see
#       the self-test docstring) show up as visible patchy/blocky artifacts
#       *in the grayscale itself*.
#   v2: kept the equal-attenuation display, and instead drew a boundary line
#       by tracking a per-ray boolean ("was liquid/gas present anywhere along
#       this ray?"), taking the deepest True point per z-column, with a
#       morphological-opening cleanup pass to drop isolated speckle. Still
#       not good enough: a boolean presence flag is maximally sensitive to
#       even a single stray sub-ray/degenerate-triangle hit, so the cleanup
#       was fighting the method's own fragility rather than fixing it.
#   v3 (current): trace the boundary from a *continuous attenuation value*
#       instead of a boolean. The ray's theoretical maximum possible
#       mu_integral (100% one phase, over the full X_WIDTH path) is a known
#       constant. Any of the *other* phases along the ray strictly lowers
#       mu_integral below that ceiling, so "is the ray entirely phase X" is
#       a continuous, averaged test, not a boolean OR of noisy samples -- one
#       stray triangle only nudges the average a little instead of flipping
#       a hard yes/no. This needs its own auxiliary attenuation computed
#       with genuinely distinct MU_SOLID_AUX/MU_LIQUID_AUX (see below) purely
#       to locate the boundary; the *displayed* image still uses equal
#       MU_SOLID/MU_LIQUID throughout.
#
# Run via pvpython (needs paraview.simple to read the OpenFOAM case):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/lateral_xray.py <case.foam> <time> <output.png>
#
from paraview.simple import *
from paraview import servermanager
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import sys
import os
import re
import time

_t_start = time.time()


def log(msg):
    print(f"[{time.time() - _t_start:7.1f}s] {msg}", flush=True)


# ── Parameters ───────────────────────────────────────────────────
FIELD_GM      = 'alpha_smoothed'  # gas/metal interface source. Recorded field
                                    # (one pass of fvc::average baked in by the
                                    # solver, see updateProps.H) -- NOT the raw
                                    # alpha.metal, and NOT re-smoothed here.
FIELD_LS      = 'epsilon1'        # liquid/solid interface source: recorded
                                    # liquid-fraction field (0=solid, 1=liquid),
                                    # already correctly clamped/blended -- no
                                    # need to recompute from T/TSolidus/TLiquidus.
ISO_THRESHOLD = 0.5               # isosurface value for both fields.
MU_GAS        = 0.0               # attenuation coeff, gas phase (per mm)
MU_SOLID      = 5.0               # attenuation coeff, solid metal (per mm) --
                                    # matches the shared MU_METAL value from
                                    # lateral_xray.py, so solid regions look
                                    # the same as in that script.
MU_LIQUID     = 5.0               # attenuation coeff, liquid metal (per mm) --
                                    # equal to MU_SOLID (see header): the melt
                                    # pool boundary is shown via the overlay
                                    # line instead of a grayscale difference.
MU_SOLID_AUX  = 5.0                # solid/liquid attenuation used *only* to
MU_LIQUID_AUX = 3.5                # locate the melt-pool boundary line (see
                                    # header, v3) -- deliberately distinct so
                                    # "is this ray 100% solid" is a meaningful,
                                    # continuous test. Never used for the
                                    # displayed image.
Y_DEPTH_MIN = 0.05e-3              # crop window: depth from nominal surface (m)
Y_DEPTH_MAX = 0.45e-3              # matches the 5um-refined mesh band -- see
                                    # tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al/
                                    # system/topoSetDict
MELT_FRONT_OFFSET = 0.05e-3        # z-window forward edge = final laser z +
                                    # this -- matches lateral_screenshot.py /
                                    # top_screenshot.py's fixed scan-length
                                    # window, so all three stacked views
                                    # share the same z-extent (previously
                                    # this script used the full mesh domain
                                    # instead, which is what a wider crop
                                    # here than the other two views traced
                                    # back to).
CROP_PADDING = 0.03e-3             # margin added behind the scan's start
X_WIDTH = 0.3e-3                   # ray path length (the "lateral"/through-
                                    # thickness beam direction): the simulation
                                    # domain is 0.64mm wide, but the real
                                    # experimental sample this compares against
                                    # was only ~0.29mm -- attenuating over the
                                    # full (wider) simulation domain overstates
                                    # the path length and crushes the whole
                                    # image's contrast. Crop the ray to the
                                    # central X_WIDTH instead of the full domain.
NY, NZ = 120, 768                  # output image resolution (rays). No longer
                                    # tied to mesh cell size -- nothing is
                                    # voxelized, this is purely an image-
                                    # quality choice.
# ─────────────────────────────────────────────────────────────────

if len(sys.argv) != 4:
    print("Usage: pvpython lateral_xray.py <case.foam> <time> <output.png>")
    sys.exit(1)

foam_file, time_value, output_png = sys.argv[1], float(sys.argv[2]), sys.argv[3]


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
reader.Createcelltopointfiltereddata = 1   # need point data for Contour --
                                            # this is the reader's own native
                                            # cell->point averaging over real
                                            # mesh connectivity (volume-
                                            # weighted over each point's true
                                            # neighboring cells), fundamentally
                                            # different from resampling onto
                                            # an independent fixed-pitch grid
                                            # (which is what produced the
                                            # sawtooth artifact in the
                                            # original voxel-resample version).
reader.UpdatePipeline(time=time_value)
log("reader loaded")

merged = MergeBlocks(Input=reader)
merged.UpdatePipeline(time=time_value)
log("blocks merged")

bounds = merged.GetDataInformation().GetBounds()
xmin, xmax, ymin, ymax, zmin, zmax = bounds
log(f"Domain bounds: x=[{xmin},{xmax}] y=[{ymin},{ymax}] z=[{zmin},{zmax}]")

_laser_table = _load_laser_time_vs_position(os.path.dirname(foam_file))
_scan_start_z = _laser_table[0][3]
_scan_end_z = _laser_table[-1][3]
zmin_crop = _scan_start_z - CROP_PADDING
zmax_crop = _scan_end_z + MELT_FRONT_OFFSET + CROP_PADDING
log(f"Crop bounds: x=[{xmin},{xmax}] (full, uncropped) "
    f"y=[{Y_DEPTH_MIN},{Y_DEPTH_MAX}] z=[{zmin_crop},{zmax_crop}]")

# ── Build the gas/metal and liquid/solid isosurfaces on the native mesh ──
gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

metal_only = Threshold(Input=merged)
metal_only.Scalars = ['POINTS', FIELD_GM]
metal_only.ThresholdRange = [ISO_THRESHOLD, 1.0]
metal_only.UpdatePipeline(time=time_value)

ls_contour = Contour(Input=metal_only)
ls_contour.ContourBy = ['POINTS', FIELD_LS]
ls_contour.Isosurfaces = [ISO_THRESHOLD]
ls_contour.UpdatePipeline(time=time_value)
ls_poly = servermanager.Fetch(ls_contour)
log(f"Liquid/solid surface: {ls_poly.GetNumberOfCells()} cells")

# NOTE: vtkModifiedBSPTree was tried here for speed (it's usually cheaper
# than vtkOBBTree for many-ray queries) but its IntersectWithLine(p1, p2,
# points, cellIds) -- the "give me every crossing" overload we need -- is
# unimplemented in this VTK build: it logs "does not yet support this
# IntersectWithLine interface" and silently returns zero intersections every
# time. vtkStaticCellLocator's *multi*-intersection overload has the exact
# same failure (verified with a standalone probe -- only vtkOBBTree
# implements it in this VTK build, but it's ~44ms/ray, over an hour/image
# for a surface this size). Its *single*-hit overload
# (IntersectWithLine(p1, p2, tol, t, x, pcoords, subId)) does work correctly,
# and vtkStaticCellLocator is built specifically for fast repeated queries
# against one static dataset -- so "all crossings" is built here by
# repeatedly asking for the next hit and nudging past it, rather than
# relying on the broken convenience method. Verified below with an explicit
# self-test before trusting it for the full run.

_NUDGE = 1e-8   # in metres -- ~500x smaller than the finest 5um mesh cell,
                # so it steps past a found hit without skipping a genuinely
                # separate adjacent crossing.
_HIT_TOL = 1e-9
_MAX_HITS_PER_RAY = 64  # generous safety cap against a pathological loop


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
    degenerate for an X-direction ray by construction, not by locator bug --
    e.g. a triangle whose 3 vertices all share the same Y sits exactly in a
    Y-constant plane, parallel to our ray direction, which is a mathematically
    degenerate "ray lies in the triangle's plane" case that any ray-triangle
    intersection correctly reports as no-hit. This happens for ls_poly in
    particular: it's contoured from a Threshold output whose cut boundary
    sits exactly on the Cartesian mesh's axis-aligned faces, and epsilon1
    can be ~0.5 right there, producing a sliver triangle lying flat in that
    cut plane. One such unlucky pick isn't a sign the locator is broken;
    only failing on *every* candidate is. (These slivers can still leave a
    handful of isolated pixels with a slightly-wrong liquid/solid split --
    acceptable here since we're integrating a continuous-tone image, not
    drawing a single fragile boundary line that would visibly kink at every
    bad triangle.)
    """
    ncells = poly.GetNumberOfCells()
    n_try = min(n_candidates, ncells)
    # Spread across the full index range, not sequential: degenerate
    # boundary-artifact triangles (see docstring) tend to be one contiguous
    # run of cell indices, so the first N cells in sequence can all be bad.
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
        f"cells -- this locator class is silently failing (same failure "
        f"mode hit earlier with vtkModifiedBSPTree/vtkStaticCellLocator's "
        f"multi-hit overload). Aborting instead of producing a wrong image."
    )


_crossing_totals = {'GM': 0, 'LS': 0, '_TEST': 0}

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

# ── Point locator on the native mesh, to seed each ray's start-of-line phase. ──
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
    liquid = metal and (ls_start_arr[pid] >= ISO_THRESHOLD)
    return metal, liquid


ys = np.linspace(Y_DEPTH_MIN, Y_DEPTH_MAX, NY)
zs = np.linspace(zmin_crop, zmax_crop, NZ)

mu_integral_avg = np.empty((NZ, NY))       # main (equal-mu) -- for display
mu_integral_aux_avg = np.empty((NZ, NY))   # auxiliary (distinct-mu) -- boundary only

# rays span only the central X_WIDTH of the domain (the "lateral"/through-
# thickness beam direction), not the full native x-extent -- see X_WIDTH above.
x_center = (xmin + xmax) / 2.0
x0, x1 = x_center - X_WIDTH / 2.0, x_center + X_WIDTH / 2.0
log(f"Ray x-extent (cropped to X_WIDTH): x=[{x0},{x1}]")

# Sub-pixel supersampling: a single ray at the exact pixel center is hostage
# to exactly where a powder particle/spatter droplet happens to sit relative
# to the fixed x0/x1 window -- a particle straddling the window edge gets an
# all-or-nothing sliver counted depending on sub-micron position, producing
# noise that's an artifact of the sampling grid, not real structure. Casting
# a small jittered grid of sub-rays per output pixel and averaging converges
# toward the true local-average attenuation instead, the way a real detector
# pixel (which integrates over a finite area) would.
N_SUB = 2  # NxN sub-ray grid per output pixel
_SUB_FRACS = [(-0.25 + 0.5 * k / (N_SUB - 1)) if N_SUB > 1 else 0.0 for k in range(N_SUB)]
dy_pitch = (Y_DEPTH_MAX - Y_DEPTH_MIN) / (NY - 1)
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
        "locator is silently failing (as vtkModifiedBSPTree did here) or "
        "the surfaces/bounds are wrong. Refusing to save a bogus image."
    )

transmission = np.exp(-mu_integral_avg)
image = transmission.T          # (y_crop, z_crop)

# Melt pool boundary line, from a continuous-attenuation-ceiling test (see
# v3 in the header) rather than a boolean presence flag: needs the
# *auxiliary* mu_integral_aux (distinct MU_SOLID_AUX/MU_LIQUID_AUX) -- its
# ceiling (100% solid, no liquid, no gas) is MU_SOLID_AUX*X_WIDTH. Any
# liquid *or* gas lowers it below that. THRESHOLD_FRAC_MELT relaxes the
# ceiling to 98% -- a ray needs at least ~2% of its path length to be
# liquid/gas before the line marks it (also windowed to near the laser --
# see below -- so this only needs to reject fine mesh-sliver noise, not
# filter distant spurious detections).
THRESHOLD_FRAC_MELT = 0.98
_X_WIDTH_MM = X_WIDTH * 1e3
_melt_ceiling = MU_SOLID_AUX * _X_WIDTH_MM

# how far behind the melt pool's leading (front) edge to keep the blue line.
# The front edge itself is derived from the laser's own known position, not
# from the (very noisy at 100%) melt detection -- see the DEBUG dump this
# replaced: at the strict 100% ceiling, mesh-sliver noise registers as
# "not pure solid" almost everywhere in the domain, all the way out to the
# z boundary, even far ahead of where the laser has actually reached. The
# laser scans in +z here (see timeVsLaserPosition), so "front" = current
# laser z position + a small forward offset (the melt pool's leading edge
# sits slightly ahead of the beam centroid), and "behind" is smaller z.
MELT_WINDOW_BEHIND_FRONT = 0.6e-3  # 0.6mm (MELT_FRONT_OFFSET, used below
                                    # too, is defined up with the z-window
                                    # params -- same "laser z + 0.05mm"
                                    # concept in both places)


def _threshold_bottom(mu_avg, ceiling, y_values, frac):
    """Per z-column (row of mu_avg), deepest y where mu_avg < frac * ceiling.

    mu_avg is (NZ, NY). NaN where the column never falls below that mark
    (i.e. that whole column is essentially at max attenuation -- undisturbed).
    """
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

_laser_z = _laser_z_at(_laser_table, time_value)  # _laser_table loaded earlier, for the z-window

_melt_front_z = _laser_z + MELT_FRONT_OFFSET
_melt_window_min_z = _melt_front_z - MELT_WINDOW_BEHIND_FRONT
liquid_bottom[(zs < _melt_window_min_z) | (zs > _melt_front_z)] = np.nan
log(f"Laser z={_laser_z*1e3:.3f}mm; melt front (laser+{MELT_FRONT_OFFSET*1e3:.2f}mm)="
    f"{_melt_front_z*1e3:.3f}mm; blue line kept for "
    f"z=[{_melt_window_min_z*1e3:.3f}, {_melt_front_z*1e3:.3f}] mm")

# Beer-Lambert transmission is long-tailed here (a handful of gas/keyhole
# pixels near 1.0 can dominate a fixed [0,1] scale). One simple rule for the
# whole image: percentile-stretch over all of it. With the ray now cropped
# to the physically-correct X_WIDTH (see above) instead of the full,
# too-wide simulation domain, the bulk floor sits at a much more moderate
# transmission (not crushed near 0), so this natural, single stretch shows
# both above-surface (plume) and below-surface (bulk/depression/liquid)
# structure without needing to special-case any region.
vmin, vmax = np.percentile(image, [1, 99])
log(f"Contrast stretch: vmin={vmin:.4f} vmax={vmax:.4f} (raw range {image.min():.4f}-{image.max():.4f})")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# True z:y physical aspect ratio (1mm of z occupies the same on-image
# distance as 1mm of y), matching how the lateral_screenshot.py ParaView
# renders are scaled -- aspect='auto' (previous behavior) stretched the
# image to fill a fixed figure box regardless of the actual z/y ratio,
# smooshing the image sideways since z (the scan direction, ~3mm here) is
# much longer than y (the depth crop, ~0.4mm). Figure size is computed from
# that same ratio so there's no letterboxing either.
DPI = 150
FIG_HEIGHT_IN = 500 / DPI
z_range_mm = (zmax_crop - zmin_crop) * 1e3
y_range_mm = (Y_DEPTH_MAX - Y_DEPTH_MIN) * 1e3
fig_width_in = FIG_HEIGHT_IN * (z_range_mm / y_range_mm)

fig, ax = plt.subplots(figsize=(fig_width_in, FIG_HEIGHT_IN), dpi=DPI)
extent = [zmin_crop * 1e3, zmax_crop * 1e3, Y_DEPTH_MAX * 1e3, Y_DEPTH_MIN * 1e3]
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
                                                                    # depth reference,
                                                                    # in place of y-axis
                                                                    # text/ticks (removed below)

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_major_locator(MultipleLocator(0.5))  # round mm values only
ax.tick_params(axis='x', labelsize=21)  # 14 * 1.5
ax.tick_params(axis='y', left=False, labelleft=False)
ax.set_xlabel('z - coord (mm)', fontsize=14)
fig.tight_layout()
fig.savefig(output_png, bbox_inches='tight')
log(f"Saved: {output_png}")
