#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/mzomoro1/bin/lbf3"
BASE_CASE="$ROOT/tutorials/laserbeamFoam/testrun3"
IMAGE="${IMAGE:-lbf3}"
END_TIME="${1:-0.00030}"
MEMORY="${2:-64g}"
RANKS=(32 48 64)

run_case() {
    local np="$1"
    local case_dir="$ROOT/tutorials/laserbeamFoam/testrun3_np${np}"

    rm -rf "$case_dir"
    cp -a "$BASE_CASE" "$case_dir"

    # Start each benchmark from identical initial conditions.
    rm -rf "$case_dir/0" "$case_dir/VTKs"
    rm -rf "$case_dir"/processor*
    rm -f "$case_dir"/log.*

    sed -i -E "s/startFrom[[:space:]]+[A-Za-z]+;/startFrom       startTime;/" "$case_dir/system/controlDict"
    sed -i -E "s/endTime[[:space:]]+[0-9.eE+-]+;/endTime         ${END_TIME};/" "$case_dir/system/controlDict"

    sed -i -E "s/numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains  ${np};/" "$case_dir/system/decomposeParDict"
    sed -i -E "s/method[[:space:]]+[A-Za-z]+;/method              scotch;/" "$case_dir/system/decomposeParDict"

    echo "=== Running NP=${np}, endTime=${END_TIME} ==="
    docker run --rm \
      -v "$case_dir:/opt/laserbeamFoam/tutorials/laserbeamFoam/testrun3" \
      --cpus "$np" --memory "$MEMORY" \
      "$IMAGE" bash -lc 'cd /opt/laserbeamFoam/tutorials/laserbeamFoam/testrun3 && ./Allrun' \
      | tee "$case_dir/log.scaling"
}

summarize_case() {
    local np="$1"
    local log_file="$ROOT/tutorials/laserbeamFoam/testrun3_np${np}/log.laserbeamFoam"

    if [[ ! -f "$log_file" ]]; then
        echo "NP=${np}: missing $log_file"
        return 0
    fi

    awk -v np="$np" '
        /^Time = / {
            t=$3+0;
            if (firstT==0) firstT=t;
            lastT=t;
        }
        /ExecutionTime =/ {
            e=$3+0;
            if (firstE==0) firstE=e;
            lastE=e;
        }
        END {
            dT=lastT-firstT;
            dE=lastE-firstE;
            if (dT>0 && dE>0) {
                simUs=dT*1e6;
                rateUsPerHour=simUs/(dE/3600.0);
                secPerUs=dE/simUs;
                printf("NP=%d  dSim=%.3f us  dWall=%.1f s  simRate=%.3f us/hour  secPerUs=%.3f\n", np, simUs, dE, rateUsPerHour, secPerUs);
            } else {
                printf("NP=%d  insufficient data for summary\n", np);
            }
        }
    ' "$log_file"
}

for np in "${RANKS[@]}"; do
    run_case "$np"
done

echo
echo "=== Scaling summary ==="
for np in "${RANKS[@]}"; do
    summarize_case "$np"
done
