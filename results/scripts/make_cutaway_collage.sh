#!/bin/bash
# Build a cutaway collage from a config file -- wrapper that runs
# make_cutaway_collage.py inside the headless ParaView image.
# See results/scripts/README_cutaway_collage.md for the full how-to.
#
# Usage (from anywhere inside the repo):
#   bash results/scripts/make_cutaway_collage.sh <config.json> [--grid] [--rerender] [--jobs N]
# Example:
#   bash results/scripts/make_cutaway_collage.sh results/scripts/collage_configs/testrun69_pores.json
set -euo pipefail
[ $# -ge 1 ] || { sed -n '2,9p' "$0"; exit 1; }

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
CONFIG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
shift
case "$CONFIG" in
  "$REPO"/*) ;;
  *) echo "Config must live inside the repo ($REPO) so the container can see it"; exit 1 ;;
esac

docker run --rm --user "$(id -u):$(id -g)" -e PYTHONUNBUFFERED=1 \
  -v "$REPO":/workspace --entrypoint /opt/paraview/bin/pvpython \
  kitware/paraview:pv-v5.8.0-osmesa-py3 \
  /workspace/results/scripts/make_cutaway_collage.py "/workspace/${CONFIG#"$REPO"/}" "$@"
