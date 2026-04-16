#!/usr/bin/env python3
"""
Generate missing .vtk.series files for all VTKs directories under
tutorials/laserbeamFoam/.

Usage:
    sudo python3 results/generate_vtk_series.py

The .vtk.series format is recognised by ParaView and maps each VTK file to its
physical simulation time.  Without it, ParaView falls back to using the frame
index (0, 1, 2 ...) as the time axis and cannot auto-group files whose names
mix decimal and scientific notation (e.g. "0.00012" vs "3e-05").

Run from the lbf3 repository root, or adjust BASE_DIR below.
"""

import os
import re
import json
import sys

BASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "tutorials", "laserbeamFoam"
)

def generate_series(vtk_dir, overwrite=False):
    vtk_files = [
        f for f in os.listdir(vtk_dir)
        if f.endswith(".vtk") and not f.endswith(".series")
    ]
    if not vtk_files:
        return

    # Collect unique laser names from filenames like "rays_laser0_<time>.vtk"
    laser_names = set()
    for f in vtk_files:
        m = re.match(r'^(rays_\w+?)_(.+)\.vtk$', f)
        if m:
            laser_names.add(m.group(1))

    for laser in sorted(laser_names):
        series_path = os.path.join(vtk_dir, laser + ".vtk.series")

        if os.path.exists(series_path) and not overwrite:
            print(f"  SKIP (exists): {series_path}")
            continue

        entries = []
        prefix = laser + "_"
        for f in vtk_files:
            if f.startswith(prefix) and f.endswith(".vtk"):
                time_str = f[len(prefix):-4]
                try:
                    entries.append({"name": f, "time": float(time_str)})
                except ValueError:
                    print(f"  WARNING: cannot parse time from '{f}', skipping")

        if not entries:
            continue

        entries.sort(key=lambda x: x["time"])

        try:
            with open(series_path, "w") as fh:
                json.dump({"file-series-version": "1.0", "files": entries}, fh, indent=2)
                fh.write("\n")
            print(f"  WROTE ({len(entries):3d} entries): {series_path}")
        except PermissionError:
            print(f"  PERMISSION DENIED: {series_path} — run with sudo chmod -R u+w first")


def main():
    base = os.path.realpath(BASE_DIR)
    if not os.path.isdir(base):
        print(f"ERROR: directory not found: {base}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning: {base}\n")

    for entry in sorted(os.listdir(base)):
        vtk_dir = os.path.join(base, entry, "VTKs")
        if not os.path.isdir(vtk_dir):
            continue
        print(f"{entry}/VTKs/")
        generate_series(vtk_dir, overwrite=False)

    print("\nDone.")


if __name__ == "__main__":
    main()
