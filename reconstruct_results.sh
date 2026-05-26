#!/usr/bin/env bash
set -euo pipefail

# Reconstruct OpenFOAM decomposed case results.
# By default, reconstructs testrun19_316L and testrun20_316L in tutorials/laserbeamFoam.
#
# Examples:
#   ./reconstruct_results.sh
#   ./reconstruct_results.sh tutorials/laserbeamFoam/testrun19_316L
#   ./reconstruct_results.sh --latest tutorials/laserbeamFoam/testrun20_316L
#   ./reconstruct_results.sh --container lbf3 --container-workspace /home/openfoam/lbf3

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_NAME=""
CONTAINER_WORKSPACE=""
SOLVER="laserbeamFoam"
LATEST_ONLY=0

DEFAULT_CASES=(
    "tutorials/laserbeamFoam/testrun19_316L"
    "tutorials/laserbeamFoam/testrun20_316L"
)

usage() {
    cat <<'EOF'
Usage:
  ./reconstruct_results.sh [options] [case_dir ...]

Options:
  --container <name>           Run reconstruction inside a running container.
  --container-workspace <path> Workspace root path inside container.
                               Default: same absolute path as host workspace.
  --solver <name>              Solver log prefix for checks (default: laserbeamFoam).
  --latest                     Reconstruct only latest time.
  -h, --help                   Show this help.

Notes:
  - If no case_dir is provided, defaults are:
      tutorials/laserbeamFoam/testrun19_316L
      tutorials/laserbeamFoam/testrun20_316L
  - Requires reconstructPar to be available either on host or inside container.
EOF
}

err() {
    echo "Error: $*" >&2
    exit 1
}

info() {
    echo "[reconstruct] $*"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || err "Command not found: $1"
}

run_host_case() {
    local case_path="$1"
    local cmd="reconstructPar"

    if [[ $LATEST_ONLY -eq 1 ]]; then
        cmd+=" -latestTime"
    fi

    info "Host: reconstructing ${case_path#$ROOT_DIR/}"
    (
        cd "$case_path"
        $cmd > log.reconstruct 2>&1
    )
    info "Host: done ${case_path#$ROOT_DIR/} (log: log.reconstruct)"
}

run_container_case() {
    local case_path="$1"
    local rel_path
    local case_in_container
    local cmd="reconstructPar"

    rel_path="$(realpath --relative-to="$ROOT_DIR" "$case_path")"
    case_in_container="$CONTAINER_WORKSPACE/$rel_path"

    if [[ $LATEST_ONLY -eq 1 ]]; then
        cmd+=" -latestTime"
    fi

    info "Container[$CONTAINER_NAME]: reconstructing $rel_path"
    docker exec "$CONTAINER_NAME" bash -lc "cd '$case_in_container' && $cmd > log.reconstruct 2>&1"
    info "Container[$CONTAINER_NAME]: done $rel_path (log: log.reconstruct)"
}

CASE_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --container)
            [[ $# -ge 2 ]] || err "--container requires a value"
            CONTAINER_NAME="$2"
            shift 2
            ;;
        --container-workspace)
            [[ $# -ge 2 ]] || err "--container-workspace requires a value"
            CONTAINER_WORKSPACE="$2"
            shift 2
            ;;
        --solver)
            [[ $# -ge 2 ]] || err "--solver requires a value"
            SOLVER="$2"
            shift 2
            ;;
        --latest)
            LATEST_ONLY=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            CASE_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ ${#CASE_ARGS[@]} -eq 0 ]]; then
    CASE_ARGS=("${DEFAULT_CASES[@]}")
fi

if [[ -n "$CONTAINER_NAME" ]]; then
    require_cmd docker
    if [[ -z "$CONTAINER_WORKSPACE" ]]; then
        CONTAINER_WORKSPACE="$ROOT_DIR"
    fi
fi

TOTAL=${#CASE_ARGS[@]}
DONE=0

for case_arg in "${CASE_ARGS[@]}"; do
    if [[ "$case_arg" = /* ]]; then
        case_path="$case_arg"
    else
        case_path="$ROOT_DIR/$case_arg"
    fi

    [[ -d "$case_path" ]] || err "Case directory not found: $case_arg"

    log_file="$case_path/log.$SOLVER"
    if [[ ! -f "$log_file" ]]; then
        info "Warning: $log_file not found (continuing with reconstructPar anyway)"
    fi

    if [[ -n "$CONTAINER_NAME" ]]; then
        run_container_case "$case_path"
    else
        require_cmd reconstructPar
        run_host_case "$case_path"
    fi

    DONE=$((DONE + 1))
    info "Progress: $DONE/$TOTAL"
done

info "All reconstructions finished."
