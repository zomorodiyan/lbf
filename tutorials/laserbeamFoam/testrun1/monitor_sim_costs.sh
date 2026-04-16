#!/usr/bin/env bash
set -u

# Live monitor for laserbeamFoam runs.
# - Follows log progress (Time, ExecutionTime)
# - Tails latest solver output lines
# - Tries to summarize component cost from profiling outputs

CASE_DIR="${1:-$(pwd)}"
REFRESH_SEC="${2:-5}"
LOG_FILE="${CASE_DIR}/log.laserbeamFoam"
PROF_DIR="${CASE_DIR}/postProcessing/profiling"

if [[ ! -d "${CASE_DIR}" ]]; then
    echo "ERROR: CASE_DIR does not exist: ${CASE_DIR}"
    exit 1
fi

if ! [[ "${REFRESH_SEC}" =~ ^[0-9]+$ ]] || [[ "${REFRESH_SEC}" -lt 1 ]]; then
    echo "ERROR: REFRESH_SEC must be a positive integer."
    exit 1
fi

extract_sum_by_pattern() {
    local file="$1"
    local pattern="$2"

    # Sum the last numeric token on lines matching pattern.
    # Works with many plain-text profiling formats.
    awk -v IGNORECASE=1 -v pat="$pattern" '
        $0 ~ pat {
            for (i = NF; i >= 1; i--) {
                if ($i ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
                    sum += $i
                    break
                }
            }
        }
        END { printf "%.10g", sum + 0 }
    ' "$file"
}

get_latest_profiling_file() {
    if [[ ! -d "${PROF_DIR}" ]]; then
        return 1
    fi

    find "${PROF_DIR}" -type f \( -name "*.dat" -o -name "*.txt" -o -name "*.log" -o -name "*.xml" \) 2>/dev/null \
        | xargs -r ls -1t 2>/dev/null \
        | head -n 1
}

print_header() {
    echo "============================================================"
    echo "Case        : ${CASE_DIR}"
    echo "Log         : ${LOG_FILE}"
    echo "Refresh     : ${REFRESH_SEC}s"
    echo "Started     : $(date '+%F %T')"
    echo "Stop monitor: Ctrl+C"
    echo "============================================================"
}

print_runtime_summary() {
    if [[ ! -f "${LOG_FILE}" ]]; then
        echo "Log status  : waiting for ${LOG_FILE}"
        return
    fi

    local t_line exec_line sim_time exec_s clock_s

    t_line="$(grep '^Time =' "${LOG_FILE}" | tail -n 1)"
    exec_line="$(grep 'ExecutionTime' "${LOG_FILE}" | tail -n 1)"

    sim_time="$(awk '/^Time =/{v=$3} END{print v+0}' "${LOG_FILE}" 2>/dev/null)"
    exec_s="$(awk '/ExecutionTime/{v=$3} END{print v+0}' "${LOG_FILE}" 2>/dev/null)"
    clock_s="$(awk '/ExecutionTime/{v=$6} END{print v+0}' "${LOG_FILE}" 2>/dev/null)"

    echo "Log status  : present"
    [[ -n "${t_line}" ]] && echo "${t_line}"
    [[ -n "${exec_line}" ]] && echo "${exec_line}"

    if [[ -n "${sim_time}" && "${sim_time}" != "0" && -n "${exec_s}" && "${exec_s}" != "0" ]]; then
        awk -v sim="${sim_time}" -v ex="${exec_s}" 'BEGIN {
            ratio = sim / ex
            printf "Sim speed   : %.6e sim-sec / wall-sec\n", ratio
        }'
    fi
}

print_profiling_summary() {
    local pf ray vof flow temp mpi total

    pf="$(get_latest_profiling_file || true)"
    if [[ -z "${pf}" ]]; then
        echo "Profiling   : no profiling file found yet (check again after first write)."
        return
    fi

    echo "Profiling   : ${pf}"

    # Broad patterns to capture common naming variants.
    ray="$(extract_sum_by_pattern "${pf}" '(laser|ray|heatSource)')"
    vof="$(extract_sum_by_pattern "${pf}" '(isoAdvection|alphaEqn|alpha\.)')"
    flow="$(extract_sum_by_pattern "${pf}" '(UEqn|pEqn|pcorr|momentum|pressure)')"
    temp="$(extract_sum_by_pattern "${pf}" '(TEqn|temperature|enthalpy)')"
    mpi="$(extract_sum_by_pattern "${pf}" '(MPI|Pstream|exchange|allreduce|reduce|scatter|gather)')"

    total="$(awk -v a="${ray}" -v b="${vof}" -v c="${flow}" -v d="${temp}" -v e="${mpi}" 'BEGIN{printf "%.10g", a+b+c+d+e}')"

    echo "Approx cost summary (sum of matched numeric entries):"
    echo "  ray tracing : ${ray}"
    echo "  vof         : ${vof}"
    echo "  fluid flow  : ${flow}"
    echo "  thermal     : ${temp}"
    echo "  mpi/comms   : ${mpi}"

    if awk -v t="${total}" 'BEGIN{exit !(t>0)}'; then
        awk -v r="${ray}" -v v="${vof}" -v f="${flow}" -v te="${temp}" -v m="${mpi}" -v t="${total}" 'BEGIN {
            printf "Approx shares: ray %.1f%% | vof %.1f%% | flow %.1f%% | thermal %.1f%% | mpi %.1f%%\n",
                100*r/t, 100*v/t, 100*f/t, 100*te/t, 100*m/t
        }'
    fi

    echo "Recent profiling lines (keywords):"
    grep -Ei 'laser|ray|isoAdvection|alphaEqn|UEqn|pEqn|TEqn|MPI|Pstream|exchange|reduce' "${pf}" | tail -n 20 || true
}

print_log_tail() {
    if [[ -f "${LOG_FILE}" ]]; then
        echo "Recent log tail:"
        tail -n 25 "${LOG_FILE}"
    fi
}

print_header

while true; do
    clear
    echo "Monitor time : $(date '+%F %T')"
    echo
    print_runtime_summary
    echo
    print_profiling_summary
    echo
    print_log_tail
    sleep "${REFRESH_SEC}"
done
