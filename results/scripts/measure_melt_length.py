# Quantifies the melt pool longitudinal length L_mp per timestep for a VDEP
# case, from the reconstructed OpenFOAM output directly (no rendering).
# Companion to measure_melt_vapor_depth.py -- same framework, same T_m
# criterion, same quasi-steady time-average -- but measures a length along
# the scan (z) direction instead of a depth below the surface.
# (Originally written by Zixun; depth band / laser-relative z-window added by
# Mehrdad to keep front-rim dynamics and spatter out of L_mp.)
#
#   L_mp  -- melt pool longitudinal length: the extent, along the laser scan
#            direction (z), of the T=T_m isotherm within the metal
#            (T_m default: liquidus, see --tmelt), restricted to the search
#            box below AND to the connected piece nearest the laser (the
#            main pool). L_mp = z_max - z_min of that piece.
#   L_all -- same, but over every isotherm piece in the search box (main
#            pool + detached sub-surface liquid pockets): the maximum
#            below-surface melt extent. L_all >= L_mp always.
#
# Search box (why): taken over the full metal region, the T=T_m isotherm also
# picks up (a) the hot, sloshing rim thrown up ahead of/around the vapor
# depression and (b) molten spatter droplets flying through the gas or landed
# elsewhere on the track -- both are metal (alpha > ISO_THRESHOLD) and above
# T_m, so either one can stretch z_max/z_min far past the actual pool. Both
# live ABOVE the nominal surface, so the main filter is a depth band:
#   --depth-min-um  (default 0)    only isotherm at/below y = SURFACE_Y + this
#   --depth-max-um  (default none) ...and above y = SURFACE_Y + this
# plus an optional laser-relative z-window (scan is +z, so "ahead" = +z):
#   --z-ahead-um    (default none) drop isotherm further than this ahead of
#                                  the laser (e.g. spatter landed in front)
#   --z-behind-um   (default none) drop isotherm further than this behind it
# Leave the z-window off unless needed: the tail of the pool can be long
# (400K/500K T0 cases), and a too-tight --z-behind-um would silently truncate
# it. Any frame whose isotherm reaches a box z-edge is flagged in the CSV
# (front_clipped/tail_clipped) so a truncation can't go unnoticed.
#
# Main-pool connectivity filter (L_mp; L_all is reported alongside without it):
# checked on testrun69 (200-400us), the depth band alone barely moves L_mp --
# the front sits a steady ~55-65um ahead of the laser and the above-surface
# rim/spatter pieces mostly lie inside the pool's own z-range. What actually
# inflated L_mp there was small DETACHED T=T_m pockets 80-130um below the
# surface, up to ~250um behind the main pool (e.g. t=321us: 726um raw vs
# 464um for the main connected pool). A depth band can't remove those, and a
# z-window tight enough to would also cut genuinely long tails, so after the
# box clip only the connected isotherm piece nearest the laser spot
# (x_laser, SURFACE_Y, z_laser) is kept.
#
# Coordinate convention (same as measure_melt_vapor_depth.py/render_view.py):
# y increases with depth into the metal, nominal surface at SURFACE_Y=0.2mm,
# so depth below surface = y - SURFACE_Y (negative = above surface).
#
# Everything else is reused unchanged from measure_melt_vapor_depth.py /
# render_view.py: FIELD_GM/FIELD_LS/ISO_THRESHOLD, the metal-side clip, the
# T-on-continuous-field contour, and T_m parsed from each case's own
# constant/transportProperties (Al 6061: Tsolidus=855K/Tliquidus=925K), not
# hardcoded.
#
# Usage: prefer the results/scripts/measure_melt_length.sh wrapper.
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from render_view import (  # noqa: E402 -- must follow sys.path fix above
    FIELD_GM, FIELD_LS, ISO_THRESHOLD,
    _load_laser_time_vs_position, _laser_z_at,
)

from paraview.simple import *  # noqa: E402
from paraview import servermanager  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SURFACE_Y = 0.2e-3  # m -- nominal flat-plate surface height, see render_view.py
EDGE_TOL = 5e-6     # m -- one finest AMR cell; isotherm this close to a box
                    # z-edge counts as clipped by the window


def read_melting_point(case_dir, choice):
    """Parse Tsolidus/Tliquidus straight from this case's own
    transportProperties (identical to measure_melt_vapor_depth.py). Filters out
    the non-physical near-zero epsilon1-band values."""
    path = os.path.join(case_dir, 'constant', 'transportProperties')
    with open(path) as f:
        text = f.read()

    def vals(name):
        return [float(v) for v in re.findall(name + r'\s+([\d.eE+-]+)\s*;', text)]

    tsol = [v for v in vals('Tsolidus') if v > 300]
    tliq = [v for v in vals('Tliquidus') if v > 300]
    if not tsol or not tliq:
        raise RuntimeError(f'Could not find a physical Tsolidus/Tliquidus (>300K) in {path}')
    Tsolidus, Tliquidus = tsol[0], tliq[0]

    if choice == 'liquidus':
        Tm = Tliquidus
    elif choice == 'solidus':
        Tm = Tsolidus
    elif choice == 'average':
        Tm = 0.5 * (Tsolidus + Tliquidus)
    else:
        Tm = float(choice)
    return Tm, Tsolidus, Tliquidus


def measure_at_time(reader, laser_table, t, Tm, box):
    """Returns (L_mp um, front um ahead of laser, tail um behind laser,
    L_all um, laser_z m, front_clipped, tail_clipped). L_mp/front/tail are
    for the connected main pool, L_all for every isotherm piece in the box.
    box = (depth_min, depth_max, z_ahead, z_behind) in metres, None =
    unbounded."""
    depth_min, depth_max, z_ahead, z_behind = box
    laser_x = laser_table[0][1]  # scan is straight along z at fixed x
    reader.UpdatePipeline(time=t)

    merged = MergeBlocks(Input=reader)
    merged.UpdatePipeline(time=t)
    xmin, xmax, ymin, ymax, zmin, zmax = merged.GetDataInformation().GetBounds()

    laser_z = _laser_z_at(laser_table, t)
    y0 = max(ymin, SURFACE_Y + depth_min) if depth_min is not None else ymin
    y1 = min(ymax, SURFACE_Y + depth_max) if depth_max is not None else ymax
    z0 = max(zmin, laser_z - z_behind) if z_behind is not None else zmin
    z1 = min(zmax, laser_z + z_ahead) if z_ahead is not None else zmax

    # search box (depth band x laser-relative z-window), full x -- same
    # Box-clip pattern as measure_melt_vapor_depth.py's z-window
    window = Clip(Input=merged)
    window.ClipType = 'Box'
    window.ClipType.Position = [xmin, y0, z0]
    window.ClipType.Length = [xmax - xmin, y1 - y0, z1 - z0]
    window.Invert = 1
    window.UpdatePipeline(time=t)

    # metal side only (same alpha criterion as h_mp)
    metal_clip = Clip(Input=window)
    metal_clip.ClipType = None
    metal_clip.Scalars = ['POINTS', FIELD_GM]
    metal_clip.Value = ISO_THRESHOLD
    metal_clip.Invert = 0
    metal_clip.UpdatePipeline(time=t)

    # T = T_m isotherm within the metal; its z-span is the melt pool length
    ls_curve = Contour(Input=metal_clip)
    ls_curve.ContourBy = ['POINTS', FIELD_LS]
    ls_curve.Isosurfaces = [Tm]
    ls_curve.UpdatePipeline(time=t)

    # every isotherm piece in the box -> L_all (max below-surface extent)
    b_all = ls_curve.GetDataInformation().GetBounds()

    # only the isotherm piece nearest the laser spot (the main pool) -> L_mp;
    # drops detached sub-surface liquid pockets behind the pool that the
    # depth band can't reach (see module header)
    pool = Connectivity(Input=ls_curve)
    pool.ExtractionMode = 'Extract Closest Point Region'
    pool.ClosestPoint = [laser_x, SURFACE_Y, laser_z]
    pool.UpdatePipeline(time=t)
    b = pool.GetDataInformation().GetBounds()

    Delete(pool)
    Delete(ls_curve)
    Delete(metal_clip)
    Delete(window)
    Delete(merged)

    L_all = (b_all[5] - b_all[4]) * 1e6 if b_all[5] >= b_all[4] else float('nan')
    iz0, iz1 = b[4], b[5]
    if iz1 < iz0:  # empty isotherm
        nan = float('nan')
        return nan, nan, nan, L_all, laser_z, 0, 0
    L_mp = (iz1 - iz0) * 1e6
    front = (iz1 - laser_z) * 1e6
    tail = (laser_z - iz0) * 1e6
    front_clipped = int(z_ahead is not None and iz1 >= z1 - EDGE_TOL)
    tail_clipped = int(z_behind is not None and iz0 <= z0 + EDGE_TOL)
    return L_mp, front, tail, L_all, laser_z, front_clipped, tail_clipped


def _um_or_none(v):
    return None if v is None else v * 1e-6


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('foam_file')
    ap.add_argument('output_csv')
    ap.add_argument('output_plot')
    ap.add_argument('--t-start', type=float, default=None, help='s, restrict analysis to t >= this')
    ap.add_argument('--t-end', type=float, default=None, help='s, restrict analysis to t <= this')
    ap.add_argument('--steady-frac', type=float, default=0.5,
                     help='fraction of the analyzed time range (from the end) treated as the '
                          'quasi-steady window for the time-average. Default 0.5 = latter half. '
                          'Ignored if --steady-start is given.')
    ap.add_argument('--steady-start', type=float, default=None,
                     help='s, absolute start time of the quasi-steady window (overrides '
                          '--steady-frac).')
    ap.add_argument('--tmelt', default='liquidus',
                     help="melt threshold T_m: 'liquidus' (default), 'solidus', 'average', "
                          "or an explicit Kelvin value -- keep consistent with h_mp/d_vd")
    ap.add_argument('--depth-min-um', type=float, default=0.0,
                     help='only count isotherm at least this far below the nominal surface '
                          '(um; negative = allow above it). Default 0 = drop everything above '
                          'the surface (front rim, spatter).')
    ap.add_argument('--depth-max-um', type=float, default=None,
                     help='only count isotherm at most this far below the nominal surface (um). '
                          'Default: no limit.')
    ap.add_argument('--z-ahead-um', type=float, default=None,
                     help='only count isotherm at most this far ahead (+z) of the laser (um). '
                          'Default: no limit.')
    ap.add_argument('--z-behind-um', type=float, default=None,
                     help='only count isotherm at most this far behind (-z) the laser (um). '
                          'Default: no limit (a tight value truncates long tails -- check the '
                          'tail_clipped column).')
    args = ap.parse_args()

    case_dir = os.path.dirname(os.path.abspath(args.foam_file))
    Tm, Tsolidus, Tliquidus = read_melting_point(case_dir, args.tmelt)
    print(f"Case: {case_dir}")
    print(f"T_m = {Tm:.1f} K  (Tsolidus={Tsolidus:.1f}K, Tliquidus={Tliquidus:.1f}K, --tmelt={args.tmelt})")
    box = (_um_or_none(args.depth_min_um), _um_or_none(args.depth_max_um),
           _um_or_none(args.z_ahead_um), _um_or_none(args.z_behind_um))
    print(f"Search box: depth [{args.depth_min_um}, {args.depth_max_um}] um below surface, "
          f"z [-{args.z_behind_um}, +{args.z_ahead_um}] um about the laser (None = unbounded)")

    laser_table = _load_laser_time_vs_position(case_dir)

    reader = OpenFOAMReader(FileName=args.foam_file)
    reader.CellArrays = [FIELD_GM, FIELD_LS]
    reader.Createcelltopointfiltereddata = 1
    reader.UpdatePipeline()
    times = sorted(reader.TimestepValues)
    if args.t_start is not None:
        times = [t for t in times if t >= args.t_start]
    if args.t_end is not None:
        times = [t for t in times if t <= args.t_end]
    if not times:
        raise SystemExit('No reconstructed timesteps in the requested range.')
    print(f"Found {len(times)} timesteps: t=[{times[0]:.6g}, {times[-1]:.6g}]s")

    rows = []
    for i, t in enumerate(times):
        L_mp, front, tail, L_all, laser_z, fc, tc = measure_at_time(
            reader, laser_table, t, Tm, box)
        rows.append((t, L_mp, front, tail, L_all, laser_z, fc, tc))
        flag = (' FRONT-CLIPPED' if fc else '') + (' TAIL-CLIPPED' if tc else '')
        print(f"[{i + 1}/{len(times)}] t={t:.6g}s  L_mp={L_mp:8.1f}um  L_all={L_all:8.1f}um  "
              f"front={front:7.1f}um  tail={tail:7.1f}um  laser_z={laser_z * 1e3:.3f}mm{flag}")

    t0, t1 = times[0], times[-1]
    steady_t0 = args.steady_start if args.steady_start is not None else t1 - args.steady_frac * (t1 - t0)
    steady_rows = [r for r in rows if r[0] >= steady_t0]

    def mean_std(col):
        v = np.array([r[col] for r in steady_rows if not np.isnan(r[col])])
        return (float(v.mean()), float(v.std())) if len(v) else (float('nan'), float('nan'))

    L_mean, L_std = mean_std(1)
    L_all_mean, L_all_std = mean_std(4)
    n_clipped = sum(1 for r in steady_rows if r[6] or r[7])

    print(f"\nQuasi-steady window: t=[{steady_t0:.6g}, {t1:.6g}]s "
          f"({len(steady_rows)} of {len(rows)} frames)")
    print(f"  L_mp  = {L_mean:.1f} +/- {L_std:.1f} um  (connected main pool)")
    print(f"  L_all = {L_all_mean:.1f} +/- {L_all_std:.1f} um  (all melt below surface)")
    if n_clipped:
        print(f"  WARNING: {n_clipped} steady-window frame(s) reach a z-window edge -- "
              f"L_mp is truncated there; widen --z-ahead-um/--z-behind-um")

    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)) or '.', exist_ok=True)
    with open(args.output_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time_s', 'L_mp_um', 'front_ahead_of_laser_um', 'tail_behind_laser_um',
                    'L_all_um', 'laser_z_mm', 'front_clipped', 'tail_clipped',
                    'in_steady_window'])
        for t, L_mp, front, tail, L_all, laser_z, fc, tc in rows:
            w.writerow([t, L_mp, front, tail, L_all, laser_z * 1e3, fc, tc,
                        int(t >= steady_t0)])
        w.writerow([])
        w.writerow(['# T_m_K', Tm])
        w.writerow(['# Tsolidus_K', Tsolidus])
        w.writerow(['# Tliquidus_K', Tliquidus])
        w.writerow(['# depth_min_um', args.depth_min_um])
        w.writerow(['# depth_max_um', args.depth_max_um])
        w.writerow(['# z_ahead_um', args.z_ahead_um])
        w.writerow(['# z_behind_um', args.z_behind_um])
        w.writerow(['# steady_window_start_s', steady_t0])
        w.writerow(['# L_mp_mean_um', L_mean])
        w.writerow(['# L_mp_std_um', L_std])
        w.writerow(['# L_all_mean_um', L_all_mean])
        w.writerow(['# L_all_std_um', L_all_std])
    print(f"Saved: {args.output_csv}")

    ts_us = np.array([r[0] for r in rows]) * 1e6
    L_arr = np.array([r[1] for r in rows])
    L_all_arr = np.array([r[4] for r in rows])
    clipped = np.array([bool(r[6] or r[7]) for r in rows])

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
    ax.axvspan(steady_t0 * 1e6, t1 * 1e6, color='gray', alpha=0.15, label='quasi-steady window')
    ax.plot(ts_us, L_all_arr, '-o', ms=3, color='#ef6c00',
            label=f'L_all (all melt below surface), mean {L_all_mean:.0f}um')
    ax.plot(ts_us, L_arr, '-o', ms=3, color='#1b5e20',
            label=f'L_mp (connected main pool), mean {L_mean:.0f}um')
    if clipped.any():
        ax.plot(ts_us[clipped], L_arr[clipped], 'x', ms=7, color='#b71c1c',
                label='reaches z-window edge (truncated)')
    ax.axhline(L_all_mean, color='#ef6c00', ls='--', lw=1)
    ax.axhline(L_mean, color='#1b5e20', ls='--', lw=1)
    ax.set_xlabel('time (us)')
    ax.set_ylabel('melt pool length along scan (um)')
    box_desc = f'depth >= {args.depth_min_um:g}um'
    if args.depth_max_um is not None:
        box_desc += f', <= {args.depth_max_um:g}um'
    if args.z_ahead_um is not None or args.z_behind_um is not None:
        box_desc += f'; z in [-{args.z_behind_um}, +{args.z_ahead_um}]um of laser'
    ax.set_title(f'{os.path.basename(case_dir)}  ({box_desc})', fontsize=9)
    ax.legend(loc='best', fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output_plot)
    print(f"Saved: {args.output_plot}")


if __name__ == '__main__':
    main()
