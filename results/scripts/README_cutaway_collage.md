# Cutaway collage: highlighting regions (pores, contact, pinch-off) in VDEP runs

`make_cutaway_collage.sh` builds a grid of side-view cutaway panels from one
simulation, one panel per timestep. In each panel you can mark **regions**:
the part of the gas/metal outline inside a box you choose is redrawn thicker
in a color you pick, and a label in the same color names it (e.g. "Contact
occurs", "Pore₁", "Pore₂"). You describe the figure in a small config file;
you don't need to edit any code.

Example output (testrun69, built from
[`collage_configs/testrun69_pores.json`](collage_configs/testrun69_pores.json)):
protrusion → contact → approaching pinch-off → pinch-off with two pores.

## What each panel shows

- The case is cut along the scan direction at `x = x_plane_um` (default −35 µm), and you look at it from the side.
  Behind the cut, the gas/metal surface is colored by x (blue = far, red = at the cut).
- On the cut face: solid metal is light gray, the mushy zone mid gray, and liquid dark gray.
  Black curves are the solidus and liquidus, read from the case's own
  `constant/transportProperties` (metal block), so other alloys work as well.
- Optional overlays: surface flow arrows (`velocity`), laser rays (`rays`), and
  flow arrows on the cut face (`section_velocity`).
- The dashed line is the nominal flat surface.

## Prerequisites

| Requirement | Notes |
|---|---|
| Docker + the ParaView image | `docker pull kitware/paraview:pv-v5.8.0-osmesa-py3`. This is **not** the `lbf3` simulation image. |
| A **reconstructed** case | Run `reconstructParMesh` and then `reconstructPar` first (see `TESTRUNS.md`). Only times that have been written to disk can be used as panels. |
| VDEP case layout | Nominal surface at y = 0.2 mm with y increasing into the metal; scan along +z; fields `alpha_smoothed` and `T`; laser path in `constant/timeVsLaserPosition`. For `rays`, you also need `VTKs/rays_laser0.vtk.series`. Any case cloned from the `tutorials/laserbeamFoam/vdep/` cases qualifies. |
| Fonts (optional) | The figures use Liberation Sans, a free substitute for Arial. The fonts are not in git. Without them you get matplotlib's default font. To install them: `sudo apt install fonts-liberation && mkdir -p results/fonts && cp /usr/share/fonts/truetype/liberation/LiberationSans-*.ttf results/fonts/` |

## Quick start

From the repo root, in WSL/Linux:

```bash
bash results/scripts/make_cutaway_collage.sh results/scripts/collage_configs/testrun69_pores.json
```

This writes `results/testrun69_cutaway_collage/testrun69_pores.png`. It takes
under a minute the first time (the 4 panels render in parallel) and a few
seconds after that, because the renders are reused (see "Re-running" below).

## Making a collage for your own simulation

**1. Copy the example config** and point it at your case:

```bash
cp results/scripts/collage_configs/testrun69_pores.json results/scripts/collage_configs/<yourcase>.json
```

Then edit `case`, `output_dir` and `name`, and list the panel `time`s you want.
To find which times exist, run `ls <case>`. Use the time directory names; if a
time isn't on disk, the nearest one is used and a note is printed. For the
first pass, delete every `regions`, `excludes` and `texts` entry.

If your melt track sits somewhere else in z, set `render.z_window_mm` (what
gets rendered) and `crop_z_mm` (what each panel shows). Both are in mm.

**2. Render with a coordinate grid** to find where your features are:

```bash
bash results/scripts/make_cutaway_collage.sh results/scripts/collage_configs/<yourcase>.json --grid
```

Open `<output_dir>/raw/*_grid.png`. Each one has a grid every 0.1 mm in z
(horizontal) and every 50 µm in y (vertical, measured from the nominal
surface, positive = deeper). For each feature you want to mark, read off a box
`[z0, z1, y0, y1]` that encloses the **outline segment** you want colored. The box
only recolors the gas/metal outline curve inside it; it doesn't draw a
rectangle.

**3. Add the regions** to each panel and run again (without `--grid`):

```json
"regions": [
  {"box": [1.72, 1.84, 60, 135], "color": "#b7f58e", "label": "Pore$_1$", "label_at": [0.13, 0.06]}
]
```

Check the log for each panel (`<output_dir>/raw/*.log`). It has one
`Highlight: ... N pts highlighted` line per box. **0 pts means the box missed
the outline**, so widen or move it. Then repeat step 3 until the figure looks
right.

**4. Remove stray fragments, if any.** If a box also catches a small
neighboring piece of outline you didn't mean to color, list a smaller box
around that piece under the panel's `excludes`. It goes back to the normal
black outline.

## Config reference

The config is JSON, and lines starting with `//` are treated as comments.
Units: **z in mm**, **y in µm measured from the nominal surface**
(negative = above the surface / gas, positive = into the metal). Text positions
(`at`, `label_at`) are **fractions of the panel**: `[0, 0]` is the bottom-left
corner and `[1, 1]` the top-right.

| Key | Default | Meaning |
|---|---|---|
| `case` | required | Case directory, relative to the repo root. A `.foam` marker is created if none exists. |
| `output_dir` | required | Where outputs go. Anything under `results/` is git-ignored, which is intended. |
| `name` | `collage` | Prefix for the output file names. |
| `render.x_plane_um` | −35 | Position of the cut plane (µm). |
| `render.y_window_um` | [−300, 300] | Vertical range to render. |
| `render.z_window_mm` | [1.0, 2.5] | Range along the scan direction to render. |
| `render.supersample` | 3 | Resolution multiplier. Higher means sharper outlines when the crop is tight. |
| `render.velocity` / `rays` / `section_velocity` | false | Overlays (see "What each panel shows"). |
| `render.color_by_y` | false | Color the surface by depth instead of x. |
| `crop_z_mm`, `crop_y_um` | full render window | Part of the render each panel shows. All panels use the same crop, so they line up. |
| `columns` | 2 | Number of panels per row in the grid. |
| `time_labels` | `relative` | `relative` (µs since panel 1), `absolute` (µs), or `none`. |
| `font_scale` | 1.0 | Multiplier for all text sizes. |
| `surface_line` | true | Draw the dashed nominal-surface line. |
| `panels[].time` | required | Simulation time (s). |
| `panels[].regions[]` | — | `box` `[z0,z1,y0,y1]`, `color` (`"#rrggbb"`), optional `label`, plus `label_at`, `label_ha` (`left`/`center`/`right`) and `label_va` (`bottom`/`center`/`top`). Without `label_at`, the label goes just below the box. |
| `panels[].excludes[]` | — | `[z0,z1,y0,y1]` boxes whose outline stays plain black. They are applied before the regions. |
| `panels[].texts[]` | — | Free text: `text`, `at`, `color`, `ha`, `va`, `size` (relative to the main font size). |
| `panels[].colorbar` | false | Place the x-color colorbar in this panel's top-right corner. |
| `panels[].scale_bar_um` | — | Draw a scale bar of this length (µm) in the top-right corner. |

Labels can use matplotlib mathtext, e.g. `"Pore$_1$"` for a subscript.

## Outputs

| File | What it is |
|---|---|
| `<output_dir>/<name>.png` | The collage. |
| `<output_dir>/<name>_panelN.png` | Each finished panel on its own, e.g. for slides. |
| `<output_dir>/raw/*_t<time>.png` | Uncropped renders. Each has an `.extent.json` (the physical coordinates the image covers), a `.flags` file (the render settings used), and a `.log`. |
| `<output_dir>/raw/*_grid.png` | Coordinate-grid copies (only with `--grid`). |
| `<output_dir>/raw/*_colorbar.png` | The colorbar image used by `colorbar: true`. |

## Re-running

- A panel is only re-rendered when its render settings change: anything under
  `render`, or that panel's `regions` boxes and colors, `excludes`, or `time`.
  Changes to labels, texts, crop, colorbar or scale bar just rebuild from the
  cached renders, which takes seconds.
- `--rerender` forces every panel to render again, e.g. after the case data changed.
- `--jobs N` sets how many panels render in parallel (default 4).

## Files

| File | Role |
|---|---|
| `make_cutaway_collage.sh` | Wrapper: runs the Python driver inside the ParaView Docker image. |
| `make_cutaway_collage.py` | Driver: renders panels, crops them, adds labels, assembles the grid. |
| `render_view.py` | Does the actual rendering (`--view=cutaway`). It also provides the other VDEP views (top/lateral/xray/transverse). |
| `collage_configs/*.json` | One config per figure. Commit yours here so others can rebuild your figure. |
