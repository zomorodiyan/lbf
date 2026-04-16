#!/bin/bash
# Benchmark script: compare 4 cores vs 24 cores
# Mesh: 20x30x40 = 24,000 cells at 5µm resolution
# Short run to endTime = 5e-7 s
set -e

. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Allow mpirun as root (Docker)
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

CASE_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$CASE_DIR"

clean_case() {
    rm -rf 0 processor* log.* [0-9]* [0-9]*e*
}

prepare_case() {
    cp -r initial 0
    blockMesh > log.blockMesh 2>&1
    setSolidFraction -compressible -alpha.phase1 "alpha.metal1" -alpha.phase2 "alpha.air" > log.setSolidFraction 2>&1
    transformPoints -rotate '((0 1 0) (0 0 1))' > log.transformPoints 2>&1
    echo "  Mesh prepared OK"
}

write_decomposeParDict() {
    local NPROCS=$1
    local N=$2
    cat > system/decomposeParDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains $NPROCS;

method          simple;

coeffs
{
    n           $N;
}
EOF
}

run_benchmark() {
    local NPROCS=$1
    local N=$2

    echo "============================================"
    echo "  Benchmark: $NPROCS cores, simple $N"
    echo "============================================"

    clean_case
    prepare_case

    write_decomposeParDict "$NPROCS" "$N"

    echo "  Decomposing..."
    decomposePar > log.decomposePar 2>&1

    echo "  Running solver on $NPROCS cores..."
    local START=$(date +%s)
    mpirun -np "$NPROCS" --oversubscribe compressibleLaserbeamFoam -parallel > "log.solver_${NPROCS}p" 2>&1
    local END=$(date +%s)
    local ELAPSED=$((END - START))

    echo "  Wall-clock time: ${ELAPSED} seconds"

    # Extract last ClockTime from log
    local CLOCK=$(grep 'ClockTime' "log.solver_${NPROCS}p" | tail -1)
    if [ -n "$CLOCK" ]; then
        echo "  Last timing line: $CLOCK"
    else
        echo "  WARNING: No ClockTime found — solver may have failed"
        echo "  Last 10 lines of log:"
        tail -10 "log.solver_${NPROCS}p"
    fi
    echo ""
}

echo ""
echo "LPBF_small_vapour benchmark — 5um mesh (24,000 cells)"
echo "endTime = 5e-7 s"
echo ""

# Run 1: 4 cores, simple (1 2 2)
run_benchmark 4 "(1 2 2)"

# Run 2: 24 cores, simple (2 3 4)
run_benchmark 24 "(2 3 4)"

echo "============================================"
echo "  Benchmark complete"
echo "============================================"
