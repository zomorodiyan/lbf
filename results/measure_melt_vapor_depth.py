# Quantifies two length scales per timestep for a VDEP case, from the
# reconstructed OpenFOAM output directly (no rendering) -- written for
# Zixun's request (2026-09-13, via Mehrdad) to measure these for the
# preheat/T0-sweep runs (testrun69/75/78/81/84):
#
#   h_mp  -- melt pool depth: nominal top surface down to the deepest point
#            of the T=T_m isotherm (T_m default: liquidus, see --tmelt).
#   d_vd  -- vapor depression depth: nominal top surface down to the lowest
#            point of the gas/metal free surface (the keyhole/depression).
#
# Both are reported as a full time series (CSV + plot) and as a time-average
# over a quasi-steady window (default: the latter half of the time range --
# override with --steady-frac/--t-start/--t-end once you've looked at the
# plot and picked a window by eye; single frames fluctuate, per Zixun's
# request to time-average).
#
# Coordinate convention (same as render_view.py, reused here rather than
# re-derived): y increases with depth into the material -- the nominal
# undisturbed top surface sits at y=SURFACE_Y=0.2mm (see render_top's own
# comment, "topoSetDict's y surface (0.2mm)"), y < SURFACE_Y is the gas/vapor
# side, y > SURFACE_Y is into the metal. So "depth below surface" for either
# quantity is just (y - SURFACE_Y). FIELD_GM='alpha_smoothed' and
# FIELD_LS='T' (not epsilon1) are also reused as-is from render_view.py --
# see that file's own header comments for why (T is the continuous field the
# solver actually solves for; epsilon1 is a derived near-step function of it
# and contours poorly).
#
# T_m is NOT hardcoded: render_view.py's own LS_TSOLIDUS/LS_TLIQUIDUS
# constants (840/867 K) are stale leftovers from an earlier AlSi10Mg-era
# template and do NOT match this repo's actual current material (Al 6061,
# Tsolidus=855K/Tliquidus=925K, confirmed straight from testrun69's own
# constant/transportProperties, per CLAUDE.md's physical-parameters section)
# -- so this script instead parses Tsolidus/Tliquidus directly out of each
# case's own transportProperties every time, and only ever uses render_view's
# LS_TSOLIDUS/LIQUIDUS-adjacent constants (FIELD_GM/FIELD_LS/ISO_THRESHOLD,
# _load_laser_time_vs_position, _laser_z_at) which don't encode a material
# temperature.
#
# Both h_mp and d_vd are measured only within a z-window centered on the
# laser's current position (--z-window-um, default 300um total width) --
# without this, a global search across the whole domain could pick up an
# unrelated feature elsewhere on the track (e.g. the trailing "protrusion"
# bump, or a pinched-off pore that happens to poke the free surface), not
# the actual melt pool/depression under the laser.
#
# Usage (run via pvpython inside the paraview image, same as render_view.py;
# prefer the results/measure_melt_vapor_depth.sh wrapper below, which
# resolves a bare case number and launches Docker for you):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/measure_melt_vapor_depth.py \
#     /workspace/<case>.foam /workspace/results/<prefix>_melt_vapor_depth.csv \
#     /workspace/results/<prefix>_melt_vapor_depth.png
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


def read_melting_point(case_dir, choice):
    """Parse Tsolidus/Tliquidus straight from this case's own
    transportProperties (not hardcoded -- see module docstring). Filters
    out non-physical near-zero values: the same file also defines a
    Tsolidus/Tliquidus pair for the epsilon1 regularization band (values
    like 1.0/10.0), unrelated to the real melting point."""
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


def _ymax_or_nan(poly_source):
    b = poly_source.GetDataInformation().GetBounds()
    ymin, ymax = b[2], b[3]
    return ymax if ymax >= ymin else float('nan')


def measure_at_time(reader, laser_table, t, Tm, z_halfwidth):
    reader.UpdatePipeline(time=t)

    merged = MergeBlocks(Input=reader)
    merged.UpdatePipeline(time=t)
    xmin, xmax, ymin, ymax, zmin, zmax = merged.GetDataInformation().GetBounds()

    laser_z = _laser_z_at(laser_table, t)
    z0 = max(zmin, laser_z - z_halfwidth)
    z1 = min(zmax, laser_z + z_halfwidth)

    window = Clip(Input=merged)
    window.ClipType = 'Box'
    window.ClipType.Position = [xmin, ymin, z0]
    window.ClipType.Length = [xmax - xmin, ymax - ymin, z1 - z0]
    window.Invert = 1
    window.UpdatePipeline(time=t)

    gm_curve = Contour(Input=window)
    gm_curve.ContourBy = ['POINTS', FIELD_GM]
    gm_curve.Isosurfaces = [ISO_THRESHOLD]
    gm_curve.UpdatePipeline(time=t)
    y_depression = _ymax_or_nan(gm_curve)

    metal_clip = Clip(Input=window)
    metal_clip.ClipType = None
    metal_clip.Scalars = ['POINTS', FIELD_GM]
    metal_clip.Value = ISO_THRESHOLD
    metal_clip.Invert = 0
    metal_clip.UpdatePipeline(time=t)

    ls_curve = Contour(Input=metal_clip)
    ls_curve.ContourBy = ['POINTS', FIELD_LS]
    ls_curve.Isosurfaces = [Tm]
    ls_curve.UpdatePipeline(time=t)
    y_melt = _ymax_or_nan(ls_curve)

    Delete(ls_curve)
    Delete(metal_clip)
    Delete(gm_curve)
    Delete(window)
    Delete(merged)

    d_vd = (y_depression - SURFACE_Y) * 1e6  # um, positive = below surface
    h_mp = (y_melt - SURFACE_Y) * 1e6
    return h_mp, d_vd, laser_z


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('foam_file')
    ap.add_argument('output_csv')
    ap.add_argument('output_plot')
    ap.add_argument('--t-start', type=float, default=None, help='s, restrict analysis to t >= this')
    ap.add_argument('--t-end', type=float, default=None, help='s, restrict analysis to t <= this')
    ap.add_argument('--steady-frac', type=float, default=0.5,
                     help='fraction of the analyzed time range (from the end) treated as the '
                          'quasi-steady window for the time-average. Default 0.5 = latter half.')
    ap.add_argument('--z-window-um', type=float, default=300.0,
                     help='full width (um) of the z-window centered on the laser position that '
                          'h_mp/d_vd are searched within, per timestep')
    ap.add_argument('--tmelt', default='liquidus',
                     help="melt threshold T_m: 'liquidus' (default), 'solidus', 'average', "
                          "or an explicit Kelvin value")
    args = ap.parse_args()

    case_dir = os.path.dirname(os.path.abspath(args.foam_file))
    Tm, Tsolidus, Tliquidus = read_melting_point(case_dir, args.tmelt)
    print(f"Case: {case_dir}")
    print(f"T_m = {Tm:.1f} K  (Tsolidus={Tsolidus:.1f}K, Tliquidus={Tliquidus:.1f}K, --tmelt={args.tmelt})")

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

    z_halfwidth = args.z_window_um * 1e-6 / 2.0
    rows = []
    for i, t in enumerate(times):
        h_mp, d_vd, laser_z = measure_at_time(reader, laser_table, t, Tm, z_halfwidth)
        rows.append((t, h_mp, d_vd, laser_z))
        print(f"[{i + 1}/{len(times)}] t={t:.6g}s  h_mp={h_mp:7.1f}um  "
              f"d_vd={d_vd:7.1f}um  laser_z={laser_z * 1e3:.3f}mm")

    t0, t1 = times[0], times[-1]
    steady_t0 = t1 - args.steady_frac * (t1 - t0)
    steady_rows = [r for r in rows if r[0] >= steady_t0]
    h_mp_vals = np.array([r[1] for r in steady_rows if not np.isnan(r[1])])
    d_vd_vals = np.array([r[2] for r in steady_rows if not np.isnan(r[2])])
    h_mp_mean = float(h_mp_vals.mean()) if len(h_mp_vals) else float('nan')
    h_mp_std = float(h_mp_vals.std()) if len(h_mp_vals) else float('nan')
    d_vd_mean = float(d_vd_vals.mean()) if len(d_vd_vals) else float('nan')
    d_vd_std = float(d_vd_vals.std()) if len(d_vd_vals) else float('nan')

    print(f"\nQuasi-steady window: t=[{steady_t0:.6g}, {t1:.6g}]s "
          f"({len(steady_rows)} of {len(rows)} frames)")
    print(f"  h_mp = {h_mp_mean:.1f} +/- {h_mp_std:.1f} um")
    print(f"  d_vd = {d_vd_mean:.1f} +/- {d_vd_std:.1f} um")

    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)) or '.', exist_ok=True)
    with open(args.output_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['time_s', 'h_mp_um', 'd_vd_um', 'laser_z_mm', 'in_steady_window'])
        for t, h_mp, d_vd, laser_z in rows:
            w.writerow([t, h_mp, d_vd, laser_z * 1e3, int(t >= steady_t0)])
        w.writerow([])
        w.writerow(['# T_m_K', Tm])
        w.writerow(['# Tsolidus_K', Tsolidus])
        w.writerow(['# Tliquidus_K', Tliquidus])
        w.writerow(['# steady_window_start_s', steady_t0])
        w.writerow(['# h_mp_mean_um', h_mp_mean])
        w.writerow(['# h_mp_std_um', h_mp_std])
        w.writerow(['# d_vd_mean_um', d_vd_mean])
        w.writerow(['# d_vd_std_um', d_vd_std])
    print(f"Saved: {args.output_csv}")

    ts_us = np.array([r[0] for r in rows]) * 1e6
    h_mp_arr = np.array([r[1] for r in rows])
    d_vd_arr = np.array([r[2] for r in rows])

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
    ax.axvspan(steady_t0 * 1e6, t1 * 1e6, color='gray', alpha=0.15, label='quasi-steady window')
    ax.plot(ts_us, h_mp_arr, '-o', ms=3, color='#b71c1c', label='h_mp (melt pool depth)')
    ax.plot(ts_us, d_vd_arr, '-o', ms=3, color='#0d47a1', label='d_vd (vapor depression depth)')
    ax.axhline(h_mp_mean, color='#b71c1c', ls='--', lw=1)
    ax.axhline(d_vd_mean, color='#0d47a1', ls='--', lw=1)
    ax.set_xlabel('time (us)')
    ax.set_ylabel('depth below nominal surface (um)')
    ax.invert_yaxis()  # deeper reads as "down" on the plot
    ax.set_title(os.path.basename(case_dir))
    ax.legend(loc='best', fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.output_plot)
    print(f"Saved: {args.output_plot}")


if __name__ == '__main__':
    main()
