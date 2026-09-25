# Aggregates the per-case measure_melt_vapor_depth.py CSVs (h_mp/d_vd time
# series + quasi-steady summary, appended as trailing "# key,value" rows)
# into one cross-case summary CSV and a T0-comparison plot -- for Zixun's
# preheat-run request (2026-09-13, via Mehrdad), covering whichever of
# testrun69/75/78/81/84 are available on this device at run time.
#
# Pure matplotlib/csv, no paraview.simple needed -- but run via pvpython
# (same kitware/paraview image as measure_melt_vapor_depth.py) since that's
# where Liberation Sans + matplotlib are already set up on this repo's
# Docker images, not the plain lbf3 image.
#
# Usage:
#   docker run --rm -v <repo>:/workspace \
#     --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/scripts/_plot_melt_vapor_depth_summary.py
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

_font_paths = [os.path.join('/workspace/results/fonts', _f) for _f in (
    'LiberationSans-Regular.ttf', 'LiberationSans-Bold.ttf',
    'LiberationSans-Italic.ttf', 'LiberationSans-BoldItalic.ttf')]
if all(os.path.exists(p) for p in _font_paths):
    fm.fontManager.ttflist.extend(fm.createFontList(_font_paths))
    plt.rcParams['font.family'] = 'Liberation Sans'

RESULTS_DIR = '/workspace/results/melt_pool_measurements'

# T0-sweep case table, per T0_sweep.md -- testrun69 is the 300K baseline
# fork (predates the T0 sweep numbering, but is the same fork template at
# the sweep's reference temperature).
CASES = [
    ('testrun75_vdep_3_Al', 100),
    ('testrun78_vdep_3_Al', 200),
    ('testrun69_vdep_3_Al', 300),
    ('testrun81_vdep_3_Al', 400),
    ('testrun84_vdep_3_Al', 500),
]


def read_summary(case):
    """Returns (kv, diff_mean_um, diff_std_um) or None if the case's CSV
    doesn't exist yet. diff = h_mp - d_vd is computed per-frame from the
    case's own steady-window rows (not from the two quantities' separately-
    reported std devs) since h_mp(t) and d_vd(t) are correlated -- the same
    frame's fluctuation moves both together -- so propagating their
    independent stds would overstate the difference's uncertainty."""
    path = os.path.join(RESULTS_DIR, f'{case}_melt_vapor_depth.csv')
    if not os.path.exists(path):
        return None
    kv = {}
    diffs = []
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        for row in reader:
            if not row:
                continue
            if row[0].startswith('#'):
                kv[row[0][2:]] = float(row[1])
                continue
            if int(row[idx['in_steady_window']]) == 1:
                h_mp = float(row[idx['h_mp_um']])
                d_vd = float(row[idx['d_vd_um']])
                if not (np.isnan(h_mp) or np.isnan(d_vd)):
                    diffs.append(h_mp - d_vd)
    diffs = np.array(diffs)
    diff_mean = float(diffs.mean()) if len(diffs) else float('nan')
    diff_std = float(diffs.std()) if len(diffs) else float('nan')
    return kv, diff_mean, diff_std


def main():
    rows = []
    for case, T0 in CASES:
        result = read_summary(case)
        if result is None:
            print(f"Skipping {case} (T0={T0}K): no {case}_melt_vapor_depth.csv found "
                  f"(not simulated/measured yet)")
            continue
        kv, diff_mean, diff_std = result
        rows.append((case, T0, kv, diff_mean, diff_std))
        print(f"{case} (T0={T0}K): h_mp={kv['h_mp_mean_um']:.1f}+/-{kv['h_mp_std_um']:.1f}um  "
              f"d_vd={kv['d_vd_mean_um']:.1f}+/-{kv['d_vd_std_um']:.1f}um  "
              f"diff(h_mp-d_vd)={diff_mean:.1f}+/-{diff_std:.1f}um")

    if not rows:
        raise SystemExit('No per-case melt_vapor_depth.csv files found -- run '
                          'measure_melt_vapor_depth.sh for at least one case first.')

    summary_csv = os.path.join(RESULTS_DIR, 'melt_vapor_depth_summary.csv')
    with open(summary_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case', 'T0_K', 'h_mp_mean_um', 'h_mp_std_um',
                     'd_vd_mean_um', 'd_vd_std_um', 'diff_mean_um', 'diff_std_um',
                     'T_m_K', 'Tsolidus_K', 'Tliquidus_K', 'steady_window_start_s'])
        for case, T0, kv, diff_mean, diff_std in rows:
            w.writerow([case, T0, kv['h_mp_mean_um'], kv['h_mp_std_um'],
                        kv['d_vd_mean_um'], kv['d_vd_std_um'], diff_mean, diff_std,
                        kv['T_m_K'], kv['Tsolidus_K'], kv['Tliquidus_K'],
                        kv['steady_window_start_s']])
    print(f"Saved: {summary_csv}")

    T0s = [T0 for _, T0, _, _, _ in rows]
    h_mp = [kv['h_mp_mean_um'] for _, _, kv, _, _ in rows]
    h_mp_err = [kv['h_mp_std_um'] for _, _, kv, _, _ in rows]
    d_vd = [kv['d_vd_mean_um'] for _, _, kv, _, _ in rows]
    d_vd_err = [kv['d_vd_std_um'] for _, _, kv, _, _ in rows]
    diff_mean = [d for _, _, _, d, _ in rows]
    diff_std = [s for _, _, _, _, s in rows]

    all_T0s = [T0 for _, T0 in CASES]
    missing_T0s = sorted(set(all_T0s) - set(T0s))
    title = 'Melt pool / vapor depression depth vs. preheat temperature (T0)'
    if missing_T0s:
        title += f'\n(T0={missing_T0s} not yet simulated)'

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    ax.errorbar(T0s, h_mp, yerr=h_mp_err, fmt='-o', capsize=3, color='#b71c1c',
                label='h_mp (melt pool depth)')
    ax.errorbar(T0s, d_vd, yerr=d_vd_err, fmt='-o', capsize=3, color='#0d47a1',
                label='d_vd (vapor depression depth)')
    ax.set_xlabel('T0, initial/reference temperature (K)')
    ax.set_ylabel('depth below nominal surface (um)')
    ax.invert_yaxis()
    ax.set_title(title)
    ax.legend(loc='best', fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plot_path = os.path.join(RESULTS_DIR, 'melt_vapor_depth_vs_T0.png')
    fig.savefig(plot_path)
    print(f"Saved: {plot_path}")

    diff_title = 'Melt-front / depression-floor separation (h_mp - d_vd) vs. T0'
    if missing_T0s:
        diff_title += f'\n(T0={missing_T0s} not yet simulated)'

    fig2, ax2 = plt.subplots(figsize=(7, 4.5), dpi=150)
    ax2.errorbar(T0s, diff_mean, yerr=diff_std, fmt='-o', capsize=3, color='#256428')
    ax2.set_xlabel('T0, initial/reference temperature (K)')
    ax2.set_ylabel('h_mp - d_vd (um)')
    ax2.set_title(diff_title)
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    diff_plot_path = os.path.join(RESULTS_DIR, 'melt_vapor_depth_diff_vs_T0.png')
    fig2.savefig(diff_plot_path)
    print(f"Saved: {diff_plot_path}")


if __name__ == '__main__':
    main()
