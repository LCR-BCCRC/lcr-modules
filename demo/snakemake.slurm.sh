#!/bin/bash

# Launches a snakefile of your choice on a Slurm cluster
# Usage: ./snakemake.slurm.sh <snakefile.smk> <target_rule> "<snakemake_flags>" [runtime_config]
# Example (conda only): ./snakemake.slurm.sh genome_Snakefile.smk all "--rerun-incomplete" runtime_config.yaml
# Example (apptainer):  ./snakemake.slurm.sh genome_Snakefile.smk all "" runtime_config.yaml
#
# snakefile.smk   The snakefile you want to run
# target_rule     The name of one of the target rules in the included Snakefiles
# snakemake_flags Additional snakemake flags (e.g. "--rerun-incomplete")
# runtime_config  Optional: path to a runtime_config.yaml file. See runtime_config.yaml
#                 for the format. Presence of an 'apptainer' section enables containers;
#                 absence keeps default conda behaviour.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

snakefile=$1
target_rule=${2:-all}
snakemake_flags=$3
runtime_config=$4

# Parse runtime config (emits empty strings when no config is provided)
eval "$(python3 "$SCRIPT_DIR/parse_runtime_config.py" "${runtime_config:-}")"

# Build conda flags (omitted in container mode — mutually exclusive)
conda_flags=()
if [ -z "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    conda_prefix="${SNAKEMAKE_CONDA_PREFIX:-$CONDA_PREFIX}"
    conda_flags+=(--use-conda --conda-prefix "$conda_prefix")
fi

# Build container flags
container_flags=()
if [ -n "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    container_flags+=("$SNAKEMAKE_CONTAINER_FLAG")
    [ -n "$SNAKEMAKE_SIF_PREFIX" ]  && container_flags+=("$SNAKEMAKE_PREFIX_FLAG" "$SNAKEMAKE_SIF_PREFIX")
    [ -n "$SNAKEMAKE_BIND_PATHS" ]  && container_flags+=("$SNAKEMAKE_ARGS_FLAG" "$SNAKEMAKE_BIND_PATHS")
fi

echo snakemake --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  "${conda_flags[@]}" "${container_flags[@]}" \
  --cluster-sync "srun -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 08:00:00" \
  -s $snakefile $target_rule

sleep 5

snakemake --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  "${conda_flags[@]}" "${container_flags[@]}" \
  --cluster-sync "srun -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 08:00:00" \
  -s $snakefile $target_rule
