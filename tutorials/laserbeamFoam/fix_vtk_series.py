"""
fix_vtk_series.py — rebuild .vtk.series files from actual files on disk.

When a laserbeamFoam parallel run is resumed, OpenFOAM overwrites the
.vtk.series file with only the frames from the new segment, discarding
entries from earlier segments.  This script scans all *.vtk files in the
VTKs/ directory, extracts timestamps from filenames, and writes a correct
series file that covers the full simulation history.

Usage:
    python3 fix_vtk_series.py [path/to/VTKs]   # default: ./VTKs
    python3 fix_vtk_series.py testrun24_316L/VTKs testrun25_316L/VTKs
"""

import json
import os
import re
import sys

# Filename pattern: {base}_{time}.vtk  e.g. rays_laser0_0.0001.vtk
_PATTERN = re.compile(r'^(.+)_([0-9eE.+\-]+)\.vtk$')


def rebuild_series(vtk_dir: str) -> None:
    vtk_dir = os.path.abspath(vtk_dir)
    if not os.path.isdir(vtk_dir):
        print(f"  [skip] not a directory: {vtk_dir}")
        return

    # Group files by base name
    groups: dict[str, list[tuple[float, str]]] = {}
    for fname in os.listdir(vtk_dir):
        m = _PATTERN.match(fname)
        if not m:
            continue
        base, time_str = m.group(1), m.group(2)
        try:
            t = float(time_str)
        except ValueError:
            continue
        groups.setdefault(base, []).append((t, fname))

    if not groups:
        print(f"  [skip] no timestamped .vtk files found in {vtk_dir}")
        return

    for base, entries in sorted(groups.items()):
        entries.sort(key=lambda x: x[0])   # sort by time
        series_path = os.path.join(vtk_dir, f"{base}.vtk.series")

        # Check what was there before
        old_count = 0
        if os.path.exists(series_path):
            try:
                with open(series_path) as fh:
                    old_count = len(json.load(fh).get("files", []))
            except Exception:
                pass

        series = {
            "file-series-version": "1.0",
            "files": [{"name": fname, "time": t} for t, fname in entries],
        }
        with open(series_path, "w") as fh:
            json.dump(series, fh, indent=2)

        print(f"  {base}.vtk.series: {old_count} -> {len(entries)} entries  [{vtk_dir}]")


def main():
    dirs = sys.argv[1:] if len(sys.argv) > 1 else ["VTKs"]
    for d in dirs:
        print(f"Rebuilding series in: {d}")
        rebuild_series(d)


if __name__ == "__main__":
    main()
