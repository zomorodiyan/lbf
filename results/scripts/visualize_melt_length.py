# Visual check of what measure_melt_length.py measures: for a few selected
# timesteps, a lateral (y-z side) view of the T=T_m isotherm, split into
#   - gray:   isotherm above the nominal surface -- discarded by the depth cut
#   - orange: isotherm pieces below the surface but NOT connected to the
#             main pool -- counted in L_all only. On testrun69 these mostly
#             sit on the floor of the trailing trench (the real free surface
#             behind the laser is ~100um below nominal), not inside solid
#   - green:  the connected main pool -- L_mp
# all projected onto the y-z plane (every vertex of the 3D isotherm surface,
# any x), plus the centerline (x = x_laser) gas/metal free surface for
# context, the laser position, and the front/tail markers + arrows for L_mp
# and L_all. Uses the same pipeline as measure_melt_length.measure_at_time
# (depth cut at SURFACE_Y, metal clip, T=T_m contour, closest-point
# connectivity seeded at the laser spot), so the drawn extents are the
# measured ones -- the per-panel numbers are printed for cross-checking
# against the CSV.
#
# Usage (paraview image, same as render_view.py):
#   docker run --rm -e PYTHONUNBUFFERED=1 -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/scripts/visualize_melt_length.py \
#     /workspace/<case>/case.foam /workspace/results/melt_pool_measurements/<out>.png \
#     <t1_us> <t2_us> ...
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from render_view import (  # noqa: E402 -- must follow sys.path fix above
    FIELD_GM, FIELD_LS, ISO_THRESHOLD,
    _load_laser_time_vs_position, _laser_z_at,
)
from measure_melt_length import SURFACE_Y, read_melting_point  # noqa: E402

from paraview.simple import *  # noqa: E402
from paraview import servermanager  # noqa: E402
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
plt.rcParams.update({'font.size': 15, 'axes.labelsize': 16, 'xtick.labelsize': 14, 'ytick.labelsize': 14})
from matplotlib.collections import LineCollection  # noqa: E402
import numpy as np  # noqa: E402


def _points(src, t):
    src.UpdatePipeline(time=t)
    data = servermanager.Fetch(src)
    if data is None or data.GetNumberOfPoints() == 0:
        return np.zeros((0, 3)), data
    return vtk_to_numpy(data.GetPoints().GetData()).copy(), data


def _segments(data, P):
    """Line segments (for LineCollection) of a polydata's line cells."""
    segs = []
    if data is None or data.GetNumberOfPoints() == 0:
        return segs
    lines = data.GetLines()
    ids = vtk_to_numpy(lines.GetData())
    i = 0
    while i < len(ids):
        n = ids[i]
        cell = ids[i + 1:i + 1 + n]
        for a, b in zip(cell[:-1], cell[1:]):
            segs.append([P[a], P[b]])
        i += n + 1
    return segs


def extract(reader, laser_table, t, Tm):
    laser_x = laser_table[0][1]
    laser_z = _laser_z_at(laser_table, t)

    merged = MergeBlocks(Input=reader)
    merged.UpdatePipeline(time=t)
    xmin, xmax, ymin, ymax, zmin, zmax = merged.GetDataInformation().GetBounds()

    # full metal isotherm, no depth cut (for the gray "discarded" part)
    metal_full = Clip(Input=merged)
    metal_full.ClipType = None
    metal_full.Scalars = ['POINTS', FIELD_GM]
    metal_full.Value = ISO_THRESHOLD
    metal_full.Invert = 0
    iso_full = Contour(Input=metal_full)
    iso_full.ContourBy = ['POINTS', FIELD_LS]
    iso_full.Isosurfaces = [Tm]
    P_full, _ = _points(iso_full, t)

    # same pipeline as measure_melt_length.measure_at_time (depth >= 0)
    window = Clip(Input=merged)
    window.ClipType = 'Box'
    window.ClipType.Position = [xmin, SURFACE_Y, zmin]
    window.ClipType.Length = [xmax - xmin, ymax - SURFACE_Y, zmax - zmin]
    window.Invert = 1
    metal = Clip(Input=window)
    metal.ClipType = None
    metal.Scalars = ['POINTS', FIELD_GM]
    metal.Value = ISO_THRESHOLD
    metal.Invert = 0
    iso = Contour(Input=metal)
    iso.ContourBy = ['POINTS', FIELD_LS]
    iso.Isosurfaces = [Tm]
    P_all, _ = _points(iso, t)
    pool = Connectivity(Input=iso)
    pool.ExtractionMode = 'Extract Closest Point Region'
    pool.ClosestPoint = [laser_x, SURFACE_Y, laser_z]
    P_main, _ = _points(pool, t)

    # centerline free surface (gas/metal boundary on the x = x_laser plane)
    sl = Slice(Input=merged)
    sl.SliceType = 'Plane'
    sl.SliceType.Origin = [laser_x, 0, 0]
    sl.SliceType.Normal = [1, 0, 0]
    fs = Contour(Input=sl)
    fs.ContourBy = ['POINTS', FIELD_GM]
    fs.Isosurfaces = [ISO_THRESHOLD]
    P_fs, fs_data = _points(fs, t)
    fs_segs = _segments(fs_data, P_fs)

    for o in (fs, sl, pool, iso, metal, window, iso_full, metal_full, merged):
        Delete(o)
    return laser_z, P_full, P_all, P_main, fs_segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('foam_file')
    ap.add_argument('output_png')
    ap.add_argument('times_us', nargs='+', type=float,
                    help='timesteps to show (us); the nearest reconstructed time is used')
    ap.add_argument('--tmelt', default='liquidus')
    args = ap.parse_args()

    case_dir = os.path.dirname(os.path.abspath(args.foam_file))
    Tm, _, _ = read_melting_point(case_dir, args.tmelt)
    laser_table = _load_laser_time_vs_position(case_dir)

    reader = OpenFOAMReader(FileName=args.foam_file)
    reader.CellArrays = [FIELD_GM, FIELD_LS]
    reader.Createcelltopointfiltereddata = 1
    reader.UpdatePipeline()
    avail = np.array(sorted(reader.TimestepValues))

    n = len(args.times_us)
    fig, axes = plt.subplots(n, 1, figsize=(12, 4.0 * n), dpi=150, squeeze=False)
    for ax, t_us in zip(axes[:, 0], args.times_us):
        t = float(avail[np.argmin(abs(avail - t_us * 1e-6))])
        laser_z, P_full, P_all, P_main, fs_segs = extract(reader, laser_table, t, Tm)

        def zd(P):  # (z rel. laser um, depth below surface um)
            return (P[:, 2] - laser_z) * 1e6, (P[:, 1] - SURFACE_Y) * 1e6

        above = P_full[P_full[:, 1] < SURFACE_Y]
        s = dict(s=1.2, lw=0, rasterized=True)
        if len(above):
            ax.scatter(*zd(above), color='#9e9e9e', alpha=0.5, **s,
                       label='isotherm above surface (discarded)')
        if len(P_all):
            ax.scatter(*zd(P_all), color='#ef6c00', alpha=0.6, **s,
                       label='detached pockets (L_all only)')
        if len(P_main):
            ax.scatter(*zd(P_main), color='#1b5e20', alpha=0.6, **s,
                       label='connected main pool (L_mp)')
        if fs_segs:
            segs = [[((a[2] - laser_z) * 1e6, (a[1] - SURFACE_Y) * 1e6),
                     ((b[2] - laser_z) * 1e6, (b[1] - SURFACE_Y) * 1e6)] for a, b in fs_segs]
            ax.add_collection(LineCollection(segs, colors='k', linewidths=0.8,
                                             label='free surface at centerline'))
        ax.axhline(0, color='#1565c0', ls=':', lw=1, label='nominal surface (depth cut)')
        ax.plot([0], [-8], marker='v', color='#c62828', ms=9, ls='none')

        # extents + arrows
        L_mp = L_all = float('nan')
        ylo = -150
        for P, col, name, yarr in ((P_all, '#ef6c00', 'L_all', ylo + 30),
                                   (P_main, '#1b5e20', 'L_mp', ylo + 72)):
            if not len(P):
                continue
            z0, z1 = (P[:, 2].min() - laser_z) * 1e6, (P[:, 2].max() - laser_z) * 1e6
            L = z1 - z0
            if name == 'L_mp':
                L_mp = L
            else:
                L_all = L
            for zz in (z0, z1):
                ax.axvline(zz, color=col, ls='--', lw=1)
            ax.annotate('', xy=(z0, yarr), xytext=(z1, yarr),
                        arrowprops=dict(arrowstyle='<->', color=col, lw=1.4))
            ax.text(0.5 * (z0 + z1), yarr - 5, f'{name} = {L:.0f} um',
                    color=col, ha='center', va='bottom', fontsize=15, fontweight='bold',
                    bbox=dict(fc='white', ec='none', pad=0.5, alpha=0.8))
        print(f"t={t * 1e6:.1f}us  L_mp={L_mp:.1f}um  L_all={L_all:.1f}um  "
              f"laser_z={laser_z * 1e3:.3f}mm")

        ax.set_xlim(-1150, 150)
        ax.set_ylim(200, ylo)  # depth increases downward
        ax.set_ylabel('depth below\nsurface (um)')
        ax.set_title(f't = {t * 1e6:.1f} us', fontsize=16, loc='left', fontweight='bold')
        ax.grid(alpha=0.25)
    axes[-1, 0].set_xlabel('z relative to laser (um)   [scan direction ->]')
    from matplotlib.lines import Line2D
    dot = dict(marker='o', ls='none', ms=9)
    handles = [
        Line2D([], [], color='#1b5e20', label='connected main pool (L_mp)', **dot),
        Line2D([], [], color='#ef6c00', label='detached pieces (L_all only)', **dot),
        Line2D([], [], color='#9e9e9e', label='isotherm above surface (discarded)', **dot),
        Line2D([], [], color='k', lw=0.8, label='free surface at centerline (x = x_laser)'),
        Line2D([], [], color='#1565c0', ls=':', label='nominal surface (depth cut)'),
        Line2D([], [], color='#c62828', marker='v', ls='none', ms=10, label='laser'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=3, fontsize=14,
               frameon=False, bbox_to_anchor=(0.5, 1.0))  # sits above the axes area
    fig.suptitle(f'{os.path.basename(case_dir)}: melt pool length measurement '
                 f'(T_m = {Tm:.0f} K, lateral projection, all x)', y=1.065, fontsize=17)
    fig.tight_layout()
    fig.savefig(args.output_png, bbox_inches='tight')
    print(f"Saved: {args.output_png}")


if __name__ == '__main__':
    main()
