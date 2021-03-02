#!/bin/bash

# Default to all targets
TARGETS=${@:-all}

# Determine the number of available cores for parallelization
NUM_CORES=$(grep -c '^processor' /proc/cpuinfo)
CPU_LOAD=$(uptime | sed -r 's/.*load average: ([0-9]+).*/\1/')
CPU_LOAD=$(expr $CPU_LOAD + 10)  # Add some buffer
CORES_AVAILABLE=$(expr $NUM_CORES - $CPU_LOAD)

# Make sure there are enough available cores to run Snakemake
if (( $CORES_AVAILABLE <= 0 )); then
    echo "There aren't enough cores available on this machine to run Snakemake."
    echo "    Number of cores: $NUM_CORES / Approximate CPU load: $CPU_LOAD"
    echo "The CPU load is calculated based on the last minute (with some buffer)."
    echo "Check out top/htop to see what other jobs are currently running."
    exit 1
fi
nice -n 10 snakemake --cores "${CORES_AVAILABLE}" --keep-going --latency-wait 120 --use-conda "$TARGETS"
