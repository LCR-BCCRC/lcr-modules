#!/bin/bash

# Launches a snakefile of your choice on the local machine
# Usage: ./run.sh <snakefile.smk> <target_rule> "<snakemake_flags>" [runtime_config]
# Example (conda only): ./run.sh genome_Snakefile.smk all "" runtime_config.conda.yaml
# Example (apptainer):  ./run.sh genome_Snakefile.smk all "" runtime_config.apptainer.yaml
#
# snakefile.smk   The snakefile you want to run
# target_rule     The name of one of the target rules in the included Snakefiles
# snakemake_flags Additional snakemake flags inside quotation marks (e.g. "--rerun-incomplete").
#                 Pass "" if you have no extra flags but need to specify a runtime_config.
# runtime_config  Optional: path to a runtime_config.yaml file. See runtime_config.yaml
#                 for the format. Presence of an 'apptainer' section enables containers;
#                 absence keeps default conda behaviour.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

snakefile=$1
TARGETS=${2:-all}
snakemake_flags=$3
runtime_config=$4

# Parse runtime config (emits empty strings when no config is provided)
eval "$(python3 "$SCRIPT_DIR/../utils/parse_runtime_config.py" "${runtime_config:-}")"

# Build conda flags (omitted in container mode — mutually exclusive)
conda_flags=()
if [ -z "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    conda_flags+=(--use-conda)
    [ -n "$SNAKEMAKE_CONDA_PREFIX" ] && conda_flags+=(--conda-prefix "$SNAKEMAKE_CONDA_PREFIX")
fi

# Build container flags
container_flags=()
if [ -n "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    container_flags+=("$SNAKEMAKE_CONTAINER_FLAG")
    [ -n "$SNAKEMAKE_SIF_PREFIX" ]  && container_flags+=("$SNAKEMAKE_PREFIX_FLAG" "$SNAKEMAKE_SIF_PREFIX")
    [ -n "$SNAKEMAKE_BIND_PATHS" ]  && container_flags+=("$SNAKEMAKE_ARGS_FLAG" "$SNAKEMAKE_BIND_PATHS")
fi

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

export XDG_CACHE_HOME="$PWD/.cache"

nice -n 10 snakemake --cores "${CORES_AVAILABLE}" $snakemake_flags -s $snakefile \
    --keep-going --latency-wait 120 "${conda_flags[@]}" \
    "${container_flags[@]}" $TARGETS
