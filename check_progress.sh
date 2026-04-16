#!/bin/bash

# Usage: ./check_progress.sh <case_directory>
# Example: ./check_progress.sh tutorials/laserbeamFoam/testrun3_np32_5um_plc

if [ -z "$1" ]; then
    echo "Usage: $0 <case_directory>"
    exit 1
fi

CASE_DIR="$1"
LOG_FILE="$CASE_DIR/log.laserbeamFoam"

if [ ! -f "$LOG_FILE" ]; then
    echo "Error: $LOG_FILE not found"
    exit 1
fi

# Extract latest values from log
LATEST_TIME=$(grep "^Time = " "$LOG_FILE" | tail -1 | awk '{print $3}')
LATEST_EXEC=$(grep "^ExecutionTime = " "$LOG_FILE" | tail -1 | awk '{print $3}')
LATEST_CLOCK=$(grep "^ClockTime = " "$LOG_FILE" | tail -1 | awk '{print $3}')

if [ -z "$LATEST_TIME" ] || [ -z "$LATEST_EXEC" ]; then
    echo "Error: Could not extract data from log file"
    exit 1
fi

# Read controlDict for target end time
ENDTIME=$(grep "^endTime" "$CASE_DIR/system/controlDict" | awk '{print $2}' | sed 's/;//')

echo "================================================================================"
echo "SIMULATION PROGRESS CHECK"
echo "================================================================================"
echo "Case: $CASE_DIR"
echo "Timestamp: $(date)"
echo ""
echo "SIMULATION TIME:"
echo "  Current:    $LATEST_TIME s"
echo "  Target:     $ENDTIME s"

progress=$(python3 -c "print(f'{float($LATEST_TIME)/float($ENDTIME)*100:.3f}')")
echo "  Progress:   $progress%"
echo ""
echo "ELAPSED TIME:"
echo "  Wall clock: $LATEST_CLOCK seconds = $(python3 -c "print(f'{float($LATEST_CLOCK)/3600:.1f} hours = {float($LATEST_CLOCK)/86400:.2f} days')")"
echo ""
echo "EXTRAPOLATION (assuming linear progress):"
python3 << PYTHON_EOF
import sys

sim_time = float("$LATEST_TIME")
wall_clock = float("$LATEST_CLOCK")
target_time = float("$ENDTIME")
progress = sim_time / target_time

if progress > 0:
    total_wall_clock = wall_clock / progress
    remaining = total_wall_clock - wall_clock
    remaining_days = remaining / 86400
    remaining_weeks = remaining_days / 7

    print(f"  Total wall clock needed: {total_wall_clock:.2e} s")
    print(f"                         = {total_wall_clock/3600:.1f} hours")
    print(f"                         = {total_wall_clock/86400:.1f} days")
    print(f"                         = {total_wall_clock/604800:.1f} weeks")
    print(f"")
    print(f"  Remaining (from now):    {remaining:.2e} s")
    print(f"                         = {remaining/3600:.1f} hours")
    print(f"                         = {remaining_days:.1f} days")
    print(f"                         = {remaining_weeks:.1f} weeks")
PYTHON_EOF
echo ""
echo "================================================================================"
echo "Tip: Run this script periodically to check if simulation stays on track."
echo "     If estimated remaining time increases significantly, physics may be changing."
echo "================================================================================"
