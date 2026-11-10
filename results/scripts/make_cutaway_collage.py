#!/usr/bin/env python3
# Config-driven cutaway collage: N panels of render_view.py --view=cutaway
# (one per timestep, same cut plane), each with optional colored "regions"
# (a recolored, thicker piece of the gas/metal outline inside a z/y box, plus
# a matching-color text label), assembled into a grid.
#
# Generalizes the one-off testrun69 pore/contact figure (results/scripts/
# _render_collage_panels.sh + _build_cutaway_collage_v1.py) so it can be
# pointed at any reconstructed VDEP-style case by writing a config file
# instead of editing code. Full how-to: results/scripts/README_cutaway_collage.md.
# Worked example: results/scripts/collage_configs/testrun69_pores.json.
#
# Run inside the ParaView image (prefer the make_cutaway_collage.sh wrapper,
# which builds this docker command for you):
#   docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 \
#     -v <repo>:/workspace --entrypoint /opt/paraview/bin/pvpython \
#     kitware/paraview:pv-v5.8.0-osmesa-py3 \
#     /workspace/results/scripts/make_cutaway_collage.py <config.json> [--grid] [--rerender] [--jobs N]
#
# Pipeline:
#   1. render  -- one render_view.py subprocess per panel (fresh pvpython per
#                 panel, as ParaView pipeline state leaks between renders in
#                 one process). Skipped when the raw panel already exists and
#                 was rendered with identical flags (see <panel>.flags), so
#                 editing only labels/text/layout re-assembles in seconds.
#   2. build   -- crop each raw panel to crop_z_mm/crop_y_um using the
#                 <panel>.extent.json render_view writes alongside it (exact
#                 physical extent of the image, so no camera constants are
#                 duplicated here), then draw surface line, border, labels,
#                 time stamp, scale bar, colorbar.
#   3. assemble -- grid of `columns` panels, white gaps.
import argparse
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))  # results/scripts -> repo
PVPYTHON = '/opt/paraview/bin/pvpython'

# Liberation Sans (metric-compatible Arial substitute, SIL OFL) -- same
# registration as render_view.py; matplotlib 3.1.1 in this image predates
# addfont(). Falls back to matplotlib's default font if the files are absent.
_font_paths = [os.path.join(REPO_ROOT, 'results', 'fonts', _f) for _f in (
    'LiberationSans-Regular.ttf', 'LiberationSans-Bold.ttf',
    'LiberationSans-Italic.ttf', 'LiberationSans-BoldItalic.ttf')]
if all(os.path.exists(p) for p in _font_paths):
    fm.fontManager.ttflist.extend(fm.createFontList(_font_paths))
    plt.rcParams['font.family'] = 'Liberation Sans'

# Style constants -- defaults reproduce the testrun69 collage. Font sizes
# are in points at DPI, for a REF_PANEL_H_PX-tall panel; they scale with the
# actual cropped panel height so text keeps the same relative size whatever
# crop window / supersample a config uses.
DPI = 100
REF_PANEL_H_PX = 1470             # testrun69 panel: 600um tall at supersample=3
MAIN_FONTSIZE = 58.5
SECONDARY_FONTSIZE = MAIN_FONTSIZE * 0.6   # scale-bar label
TIME_LABEL_COLOR = '#333333'
SURFACE_LINE_STYLE = dict(color='#333333', linewidth=3.0, linestyle=(0, (6, 4)))
BORDER_STYLE = dict(fill=False, edgecolor='black', linewidth=1.25, clip_on=False)
SCALE_BAR_COLOR = '#444444'
SCALE_BAR_LINEWIDTH = 16
COLORBAR_MARGIN = 0.03            # axes fraction, from the top-right corner
COLORBAR_WIDTH_FRAC = 0.512
GAP_PX = 8                        # white gap between panels in the grid
DEFAULT_REGION_COLOR = '#ffffff'
EXCLUDE_COLOR = 'black'           # render_view's "exclude" highlight color


def log(msg):
    print(msg, flush=True)


def load_config(path):
    """JSON, plus whole-line // comments (stripped before parsing) so the
    example config can explain itself."""
    with open(path) as f:
        text = ''.join(line for line in f if not line.lstrip().startswith('//'))
    return json.loads(text)


def repo_path(p):
    return p if os.path.isabs(p) else os.path.join(REPO_ROOT, p)


def case_times(case_dir):
    """Every numeric time directory in the (reconstructed) case."""
    times = []
    for name in os.listdir(case_dir):
        try:
            t = float(name)
        except ValueError:
            continue
        if t > 0 and os.path.isdir(os.path.join(case_dir, name)):
            times.append((t, name))
    return sorted(times)


def snap_time(requested, times):
    """Nearest written time directory -- the panel times in a config are
    usually typed by hand, and ParaView silently snaps anyway, so do it
    explicitly and say so."""
    if not times:
        sys.exit("No reconstructed time directories found -- run reconstructParMesh "
                 "then reconstructPar first (see TESTRUNS.md).")
    t, name = min(times, key=lambda tn: abs(tn[0] - requested))
    if abs(t - requested) > 1e-9:
        log(f"  time {requested:g} not written; using nearest {name}")
    return name


def render_flags(cfg, panel):
    r = cfg.get('render', {})
    flags = ['--view=cutaway', '--top-crop-frac=0',
             f"--supersample={r.get('supersample', 3)}",
             f"--x-plane-um={r.get('x_plane_um', -35)}"]
    y0, y1 = r.get('y_window_um', [-300, 300])
    flags += [f'--y-min-um={y0}', f'--y-max-um={y1}']
    if 'z_window_mm' in r:
        z0, z1 = r['z_window_mm']
        flags += [f'--z-min-mm={z0}', f'--z-max-mm={z1}']
    for key, flag in (('velocity', '--velocity'), ('rays', '--rays'),
                      ('section_velocity', '--section-velocity'), ('color_by_y', '--color-y'),
                      ('dotted_solidus', '--dotted-solidus')):
        if r.get(key):
            flags.append(flag)
    # Excludes first: render_view processes boxes in order and an exclude
    # only carves out of a box listed after it.
    for box in panel.get('excludes', []):
        flags.append('--highlight=' + ','.join(str(v) for v in box) + ',' + EXCLUDE_COLOR)
    for reg in panel.get('regions', []):
        flags.append('--highlight=' + ','.join(str(v) for v in reg['box']) + ','
                     + reg.get('color', DEFAULT_REGION_COLOR))
    return flags


def render_panel(foam, tname, raw_png, flags, grid, rerender):
    flags_file = os.path.splitext(raw_png)[0] + '.flags'
    extent_file = os.path.splitext(raw_png)[0] + '.extent.json'
    flags_str = ' '.join(flags)
    if (not rerender and not grid and os.path.exists(raw_png) and os.path.exists(extent_file)
            and os.path.exists(flags_file) and open(flags_file).read() == flags_str):
        log(f"  cached: {os.path.basename(raw_png)}")
        return
    cmd = [PVPYTHON, os.path.join(SCRIPT_DIR, 'render_view.py')] + flags
    if grid:
        cmd.append('--grid')
    cmd += [foam, tname, raw_png]
    log_path = os.path.splitext(raw_png)[0] + '.log'
    log(f"  rendering {os.path.basename(raw_png)} (log: {log_path})")
    with open(log_path, 'w') as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"render_view.py failed for {raw_png} -- see {log_path}")
    with open(flags_file, 'w') as f:
        f.write(flags_str)


def build_panel(cfg, panel, raw_png, ref_time_s, out_png):
    with open(os.path.splitext(raw_png)[0] + '.extent.json') as f:
        ext = json.load(f)
    img = mpimg.imread(raw_png)
    h0, w0 = img.shape[:2]
    r = cfg.get('render', {})
    crop_z = cfg.get('crop_z_mm', [ext['z_left_mm'], ext['z_right_mm']])
    crop_y = cfg.get('crop_y_um', r.get('y_window_um', [-300, 300]))

    def col(z_mm):
        return round((z_mm - ext['z_left_mm']) / (ext['z_right_mm'] - ext['z_left_mm']) * (w0 - 1))

    def row(y_um):
        return round((y_um - ext['y_top_offset_um'])
                     / (ext['y_bottom_offset_um'] - ext['y_top_offset_um']) * (h0 - 1))

    c0, c1 = max(0, col(crop_z[0])), min(w0 - 1, col(crop_z[1]))
    r0, r1 = max(0, row(crop_y[0])), min(h0 - 1, row(crop_y[1]))
    cropped = img[r0:r1 + 1, c0:c1 + 1]
    h, w = cropped.shape[:2]
    fs = MAIN_FONTSIZE * cfg.get('font_scale', 1.0) * h / REF_PANEL_H_PX

    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    # y axis points down (depth), matching the render: top = crop_y[0].
    ax.imshow(cropped, extent=(crop_z[0], crop_z[1], crop_y[1], crop_y[0]), aspect='auto')
    ax.set_xlim(crop_z)
    ax.set_ylim(crop_y[1], crop_y[0])
    if cfg.get('surface_line', True):
        ax.axhline(0, **SURFACE_LINE_STYLE)  # y offset 0 = nominal flat surface
    ax.add_patch(plt.Rectangle((crop_z[0], crop_y[0]), crop_z[1] - crop_z[0],
                               crop_y[1] - crop_y[0], **BORDER_STYLE))

    mode = cfg.get('time_labels', 'relative')
    if mode != 'none':
        t_us = round(ext['time_s'] * 1e6) - (round(ref_time_s * 1e6) if mode == 'relative' else 0)
        ax.text(0.015, 0.97, f"{t_us} μs", transform=ax.transAxes, fontsize=fs,
                fontweight='bold', color=TIME_LABEL_COLOR, ha='left', va='top')

    def put_text(text, at, color, ha='left', va='bottom', size=1.0):
        ax.text(at[0], at[1], text, transform=ax.transAxes, fontsize=fs * size,
                fontweight='bold', color=color, ha=ha, va=va)

    for reg in panel.get('regions', []):
        if not reg.get('label'):
            continue
        if 'label_at' in reg:
            put_text(reg['label'], reg['label_at'], reg.get('color', DEFAULT_REGION_COLOR),
                     reg.get('label_ha', 'left'), reg.get('label_va', 'bottom'))
        else:
            # Default: centered just below the box (y grows downward).
            z0, z1, y0, y1 = reg['box']
            ax.text((z0 + z1) / 2, y1 + 10, reg['label'], fontsize=fs, fontweight='bold',
                    color=reg.get('color', DEFAULT_REGION_COLOR), ha='center', va='top')

    for t in panel.get('texts', []):
        put_text(t['text'], t['at'], t.get('color', '#222222'), t.get('ha', 'left'),
                 t.get('va', 'bottom'), t.get('size', 1.0))

    if panel.get('scale_bar_um'):
        length_mm = panel['scale_bar_um'] * 1e-3
        z1 = crop_z[1] - 0.05 * (crop_z[1] - crop_z[0])
        z0 = z1 - length_mm
        y_bar = crop_y[0] + 0.1 * (crop_y[1] - crop_y[0])
        ax.plot([z0, z1], [y_bar, y_bar], color=SCALE_BAR_COLOR,
                linewidth=SCALE_BAR_LINEWIDTH * h / REF_PANEL_H_PX, solid_capstyle='butt')
        ax.text(z0 - 0.02 * (crop_z[1] - crop_z[0]), y_bar, f"{panel['scale_bar_um']:g} μm",
                color=SCALE_BAR_COLOR, fontsize=fs * 0.6, fontweight='bold', ha='right', va='center')

    if panel.get('colorbar'):
        # render_view saves the cutaway colorbar next to the raw panel
        # (transparent background) without pasting it -- inset it here.
        cb_img = mpimg.imread(re.sub(r'_t[\d.eE+-]+\.png$', '_colorbar.png', raw_png))
        cb_h_frac = COLORBAR_WIDTH_FRAC / (cb_img.shape[1] / cb_img.shape[0]) * (w / h)
        cb_ax = ax.inset_axes([1 - COLORBAR_MARGIN - COLORBAR_WIDTH_FRAC,
                               1 - COLORBAR_MARGIN - cb_h_frac, COLORBAR_WIDTH_FRAC, cb_h_frac])
        cb_ax.imshow(cb_img)
        cb_ax.patch.set_alpha(0)
        cb_ax.axis('off')

    ax.axis('off')
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    log(f"  built {os.path.basename(out_png)} ({w}x{h}px)")


def assemble(panel_pngs, columns, out_png):
    imgs = [mpimg.imread(p)[..., :3] for p in panel_pngs]
    h = max(i.shape[0] for i in imgs)
    w = max(i.shape[1] for i in imgs)
    rows = math.ceil(len(imgs) / columns)
    canvas = np.ones((rows * h + (rows - 1) * GAP_PX, columns * w + (columns - 1) * GAP_PX, 3),
                     dtype=imgs[0].dtype)
    for k, im in enumerate(imgs):
        rr, cc = divmod(k, columns)
        y, x = rr * (h + GAP_PX), cc * (w + GAP_PX)
        canvas[y:y + im.shape[0], x:x + im.shape[1]] = im
    mpimg.imsave(out_png, canvas)
    log(f"Saved collage: {out_png}")


def main():
    ap = argparse.ArgumentParser(
        description="Config-driven cutaway collage -- see results/scripts/README_cutaway_collage.md")
    ap.add_argument('config', help="collage config (JSON, // comment lines allowed)")
    ap.add_argument('--grid', action='store_true',
                    help="force re-render and also save <panel>_grid.png with a z(mm)/"
                         "y(um-from-surface) grid -- use it to read off region box coordinates")
    ap.add_argument('--rerender', action='store_true', help="ignore cached raw panels")
    ap.add_argument('--jobs', type=int, default=4, help="parallel renders (default 4)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    case_dir = repo_path(cfg['case'])
    foams = [f for f in os.listdir(case_dir) if f.endswith('.foam')]
    foam = os.path.join(case_dir, foams[0] if foams else 'case.foam')
    if not foams:
        open(foam, 'a').close()  # ParaView only needs the marker to exist
    out_dir = repo_path(cfg['output_dir'])
    os.makedirs(os.path.join(out_dir, 'raw'), exist_ok=True)
    name = cfg.get('name', 'collage')

    log(f"Case: {case_dir}")
    times = case_times(case_dir)
    panels = cfg['panels']
    raw_pngs = []
    tnames = []
    jobs = []
    for i, p in enumerate(panels, 1):
        tname = snap_time(p['time'], times)
        # render_view derives the colorbar's filename by stripping this
        # trailing _t<time> (-> <name>_panel<i>_colorbar.png).
        raw = os.path.join(out_dir, 'raw', f"{name}_panel{i}_t{tname}.png")
        raw_pngs.append(raw)
        tnames.append(tname)
        jobs.append((foam, tname, raw, render_flags(cfg, p), args.grid, args.rerender))

    log(f"Rendering {len(jobs)} panels ({args.jobs} at a time)...")
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for fut in [ex.submit(render_panel, *j) for j in jobs]:
            fut.result()

    ref_time_s = float(tnames[0])  # 'relative' time labels count from panel 1
    built = []
    for i, (p, raw) in enumerate(zip(panels, raw_pngs), 1):
        out = os.path.join(out_dir, f"{name}_panel{i}.png")
        build_panel(cfg, p, raw, ref_time_s, out)
        built.append(out)
    assemble(built, cfg.get('columns', 2), os.path.join(out_dir, f"{name}.png"))


if __name__ == '__main__':
    main()
