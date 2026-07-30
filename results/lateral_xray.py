# Lateral (through-thickness) X-ray-style projection of the melt pool.
#
# Mimics in-situ synchrotron X-ray radiography of LPBF: a virtual "beam"
# passes laterally through the sample thickness (x, the direction across
# the track) and the attenuation is integrated along that path for every
# point in the build(y)-scan(z) plane, producing a grayscale side view of
# the melt pool / keyhole, cropped to the finely-meshed region near the
# surface (y) and the middle of the scan track (z).
#
# v2: surface-based ray tracing (no voxel resampling). The original version
# resampled the crop region onto a uniform ResampleToImage grid and
# thresholded per voxel -- that binarizes against the mesh's real (jagged)
# cell geometry at AMR refinement-transition cells (see topoSetDict), which
# produced a sawtooth/staircase artifact no choice of grid pitch could fix.
# This version instead builds a true marching-cubes isosurface on the native
# mesh from the solver's already-recorded, pre-smoothed alpha_smoothed field,
# and computes each pixel's attenuation by exact ray/surface intersection
# instead of voxel counting.
#
# v2 profiling (first run): the gas/metal surface has 288k triangles (this
# case has discrete powder particles, so the surface includes every particle
# in the domain, not just the melt pool) -- ~44 ms/ray against it via
# vtkOBBTree, i.e. over an hour per image. Fix: vtkStaticCellLocator (built
# for fast repeated queries against one static mesh, unlike vtkOBBTree whose
# oriented-box hierarchy degrades badly against many small, scattered,
# differently-oriented surfaces like a powder bed) plus unbuffered timing
# output. See the comment above crossings() below for a correctness gotcha
# hit along the way: this locator's "give me every crossing" convenience
# overload is silently broken in this VTK build, so crossings() gets all
# crossings via repeated single-hit queries instead.
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
ISO_THRESHOLD = 0.5               # isosurface value.
MU_GAS        = 0.0               # attenuation coeff, gas phase (per mm)
MU_METAL      = 5.0               # attenuation coeff, metal phase (per mm).
                                    # Solid and liquid metal share this value:
                                    # the solver treats them as one
                                    # incompressible phase (same rho), so
                                    # there's no real density jump at melting
                                    # to justify different coefficients here.
Y_DEPTH_MIN = 0.05e-3              # crop window: depth from nominal surface (m)
Y_DEPTH_MAX = 0.45e-3              # matches the 5um-refined mesh band -- see
                                    # tutorials/laserbeamFoam/vdep/testrun64_vdep_3_Al/
                                    # system/topoSetDict
Z_WINDOW_FRAC = 1.0                # fraction of the z (scan) domain range to
                                    # keep -- 1.0 = full domain, no z crop.
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
reader.CellArrays = [FIELD_GM]
reader.Createcelltopointfiltereddata = 1   # need point data for Contour --
                                            # this is the reader's own native
                                            # cell->point averaging over real
                                            # mesh connectivity (volume-
                                            # weighted over each point's true
                                            # neighboring cells), fundamentally
                                            # different from resampling onto
                                            # an independent fixed-pitch grid
                                            # (which is what produced the
                                            # sawtooth artifact in v1).
reader.UpdatePipeline(time=time_value)
log("reader loaded")

merged = MergeBlocks(Input=reader)
merged.UpdatePipeline(time=time_value)
log("blocks merged")

bounds = merged.GetDataInformation().GetBounds()
xmin, xmax, ymin, ymax, zmin, zmax = bounds
log(f"Domain bounds: x=[{xmin},{xmax}] y=[{ymin},{ymax}] z=[{zmin},{zmax}]")

z_trim = (1.0 - Z_WINDOW_FRAC) / 2.0
zmin_crop = zmin + z_trim * (zmax - zmin)
zmax_crop = zmax - z_trim * (zmax - zmin)
log(f"Crop bounds: x=[{xmin},{xmax}] (full, uncropped) "
    f"y=[{Y_DEPTH_MIN},{Y_DEPTH_MAX}] z=[{zmin_crop},{zmax_crop}]")

# ── Build the gas/metal isosurface on the native mesh ──
# (No bounding-box clip before this: it turned out unnecessary for
# performance -- the vtkStaticCellLocator swap below made the full-domain
# surface fast enough on its own -- and a first attempt at clipping here
# produced an incomplete/flattened image, most likely an inverted Box clip
# silently discarding the region of interest. Not worth the correctness risk
# for no performance benefit, so building straight from the full merged mesh.)
gm_contour = Contour(Input=merged)
gm_contour.ContourBy = ['POINTS', FIELD_GM]
gm_contour.Isosurfaces = [ISO_THRESHOLD]
gm_contour.UpdatePipeline(time=time_value)
gm_poly = servermanager.Fetch(gm_contour)
log(f"Gas/metal surface: {gm_poly.GetNumberOfCells()} cells")

# NOTE: vtkModifiedBSPTree was tried here for speed (it's usually cheaper
# than vtkOBBTree for many-ray queries) but its IntersectWithLine(p1, p2,
# points, cellIds) -- the "give me every crossing" overload we need -- is
# unimplemented in this VTK build: it logs "does not yet support this
# IntersectWithLine interface" and silently returns zero intersections every
# time. That produced a fast-but-wrong run (every ray came out as pure gas
# or pure metal, no interior transitions -- see git history for that broken
# version). vtkStaticCellLocator's *multi*-intersection overload has the
# exact same "does not yet support this interface" failure (verified with a
# standalone probe against vtkModifiedBSPTree/vtkStaticCellLocator/
# vtkCellLocator -- only vtkOBBTree implements it in this VTK build). Its
# *single*-hit overload (IntersectWithLine(p1, p2, tol, t, x, pcoords, subId))
# does work correctly, though, and vtkStaticCellLocator is built specifically
# for fast repeated queries against one static dataset -- so "all crossings"
# is built here by repeatedly asking for the next hit and nudging past it,
# rather than relying on the broken convenience method. Verified below with
# an explicit self-test before trusting it for the full run.

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
    intersection correctly reports as no-hit. One such unlucky pick isn't a
    sign the locator is broken; only failing on *every* candidate is.
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


_crossing_totals = {'GM': 0, '_TEST': 0}

gm_tree = vtk.vtkStaticCellLocator()
gm_tree.SetDataSet(gm_poly)
gm_tree.BuildLocator()
log("gm_tree built")
_self_test_locator(gm_tree, gm_poly, "gm_tree")

# ── Point locator on the native mesh, to seed each ray's start-of-line phase. ──
merged_data = servermanager.Fetch(merged)
gm_start_arr = vtk_to_numpy(merged_data.GetPointData().GetArray(FIELD_GM))
point_locator = vtk.vtkPointLocator()
point_locator.SetDataSet(merged_data)
point_locator.BuildLocator()
log("point locator built")


def start_state(x0, y0, z0):
    pid = point_locator.FindClosestPoint((x0, y0, z0))
    return gm_start_arr[pid] >= ISO_THRESHOLD


ys = np.linspace(Y_DEPTH_MIN, Y_DEPTH_MAX, NY)
zs = np.linspace(zmin_crop, zmax_crop, NZ)

mu_integral_avg = np.empty((NZ, NY))

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
    metal = start_state(x0, y0, z0)
    xs = []
    crossings(gm_tree, (x0, y0, z0), (x1, y0, z0), 'GM', xs)
    xs.sort(key=lambda t: t[0])

    mu_integral = 0.0
    x_prev = x0
    for x_c, tag in xs + [(x1, None)]:
        seg_len_mm = (x_c - x_prev) * 1e3
        mu_integral += (MU_METAL if metal else MU_GAS) * seg_len_mm
        metal = not metal
        x_prev = x_c
    return mu_integral


log(f"starting ray loop: {NY}x{NZ} pixels x {N_SUB * N_SUB} sub-rays each")
_ray_t0 = time.time()

for zi, z0 in enumerate(zs):
    for yi, y0 in enumerate(ys):
        mu_sum = 0.0
        for dzf in _SUB_FRACS:
            zc = z0 + dzf * dz_pitch
            for dyf in _SUB_FRACS:
                yc = y0 + dyf * dy_pitch
                mu_sum += _trace_ray(yc, zc)

        mu_integral_avg[zi, yi] = mu_sum / (N_SUB * N_SUB)

    if zi == 0:
        per_row = time.time() - _ray_t0
        eta = per_row * (NZ - 1)
        log(f"row 0/{NZ} done in {per_row:.2f}s -> ETA {eta / 60:.1f} min for remaining rows")
    elif zi % 64 == 0:
        elapsed = time.time() - _ray_t0
        eta = elapsed / (zi + 1) * (NZ - zi - 1)
        log(f"row {zi}/{NZ} (ETA {eta / 60:.1f} min)")

log("ray loop done")
log(f"Total crossings found: GM={_crossing_totals['GM']} across {NY * NZ} rays")
if _crossing_totals['GM'] == 0:
    raise RuntimeError(
        "Zero gas/metal crossings found across the entire ray loop -- the "
        "locator is silently failing (as vtkModifiedBSPTree did here) or "
        "the surfaces/bounds are wrong. Refusing to save a bogus image."
    )

transmission = np.exp(-mu_integral_avg)
image = transmission.T          # (y_crop, z_crop)

# Vapor-depression boundary: the ray's theoretical maximum mu_integral (if it
# were 100% metal, no gas at all, over the full X_WIDTH path) is a known
# constant, MU_METAL*X_WIDTH. Any gas along the ray strictly lowers
# mu_integral below that ceiling -- a continuous-valued test, not a boolean
# presence flag, so a single stray sub-ray/triangle artifact only nudges the
# *averaged* mu_integral a little instead of flipping a hard yes/no. Rather
# than requiring mu_integral to sit at (numerically) 100% of the ceiling to
# count as "no gas", THRESHOLD_FRAC relaxes that to 78%: a ray needs at least
# ~22% of its path length to be gas before the line marks it, so trace
# amounts of gas (a sliver of a spatter droplet, mesh noise) don't trigger
# the line. Per z-column, trace the deepest y that still falls below that
# mark -- that's the vapor depression's bottom.
THRESHOLD_FRAC = 0.78
_gas_ceiling = MU_METAL * (X_WIDTH * 1e3)

# z-window: fixed lower bound (0.5mm) up to whichever is larger of a fixed
# 2.5mm or the laser's own current position + a small forward margin -- see
# lateral_xray_liquid_solid.py for the same logic (kept in sync here).
WHITE_WINDOW_Z_MIN = 0.5e-3            # 0.5mm
WHITE_WINDOW_Z_MAX_FIXED = 2.5e-3      # 2.5mm
WHITE_WINDOW_LASER_OFFSET = 0.1e-3     # laser z + 0.1mm


def _threshold_bottom(mu_avg, ceiling, y_values):
    """Per z-column (row of mu_avg), deepest y where mu_avg < THRESHOLD_FRAC * ceiling.

    mu_avg is (NZ, NY). NaN where the column never falls below that mark
    (i.e. that whole column is essentially at max attenuation -- no
    depression).
    """
    below = mu_avg < (THRESHOLD_FRAC * ceiling)
    bottom = np.full(mu_avg.shape[0], np.nan)
    for zi in range(mu_avg.shape[0]):
        true_idx = np.nonzero(below[zi])[0]
        if true_idx.size:
            bottom[zi] = y_values[true_idx.max()]
    return bottom


gas_bottom = _threshold_bottom(mu_integral_avg, _gas_ceiling, ys)

_laser_table = _load_laser_time_vs_position(os.path.dirname(foam_file))
_laser_z = _laser_z_at(_laser_table, time_value)
_white_window_max_z = max(WHITE_WINDOW_Z_MAX_FIXED, _laser_z + WHITE_WINDOW_LASER_OFFSET)
gas_bottom[(zs < WHITE_WINDOW_Z_MIN) | (zs > _white_window_max_z)] = np.nan
log(f"White line z window: [{WHITE_WINDOW_Z_MIN*1e3:.3f}, {_white_window_max_z*1e3:.3f}] mm")

# Beer-Lambert transmission is long-tailed here (a handful of gas/keyhole
# pixels near 1.0 can dominate a fixed [0,1] scale). One simple rule for the
# whole image: percentile-stretch over all of it. With the ray now cropped
# to the physically-correct X_WIDTH (see above) instead of the full,
# too-wide simulation domain, the bulk floor sits at a much more moderate
# transmission (not crushed near 0), so this natural, single stretch shows
# both above-surface (plume) and below-surface (bulk/depression) structure
# without needing to special-case either region.
vmin, vmax = np.percentile(image, [1, 99])
log(f"Contrast stretch: vmin={vmin:.4f} vmax={vmax:.4f} (raw range {image.min():.4f}-{image.max():.4f})")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(7, 4), dpi=250)
extent = [zmin_crop * 1e3, zmax_crop * 1e3, Y_DEPTH_MAX * 1e3, Y_DEPTH_MIN * 1e3]
ax.imshow(
    image,
    cmap='gray',
    vmin=vmin, vmax=vmax,
    extent=extent,
    aspect='auto',
    interpolation='bilinear',
)
ax.plot(zs * 1e3, gas_bottom * 1e3, linestyle=':', color='white', linewidth=1.2,
        label='vapor depression bottom (gas/metal)')
ax.set_xlabel('z - scan direction (mm)')
ax.set_ylabel('y - depth from surface (mm)')
ax.set_title(f'lateral X-ray-style projection, t={time_value*1e6:.0f} us')
fig.tight_layout()
fig.savefig(output_png)
log(f"Saved: {output_png}")
