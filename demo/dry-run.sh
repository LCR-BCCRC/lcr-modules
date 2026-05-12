#!/bin/bash

# Launches a snakefile of your choice in dry run mode (for debugging)
# Usage: ./dry-run.sh <snakefile.smk> <target_rule> "<snakemake_flags>" [runtime_config]
# Example (conda only): ./dry-run.sh genome_Snakefile.smk all "" runtime_config.yaml
# Example (apptainer):  ./dry-run.sh genome_Snakefile.smk all "" runtime_config.yaml
#
# snakefile.smk   The snakefile you want to run
# target_rule     The name of one of the target rules in the included Snakefiles
# snakemake_flags Additional snakemake flags, inside quotation marks
# runtime_config  Optional: path to a runtime_config.yaml file. Pass the same file
#                 you use in run.sh / snakemake.slurm.sh so the dry run reflects
#                 the exact command that will be executed.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

snakefile=$1
TARGETS=${2:-all}
snakemake_flags=$3
runtime_config=$4

# Parse runtime config (emits empty strings when no config is provided)
eval "$(python3 "$SCRIPT_DIR/parse_runtime_config.py" "${runtime_config:-}")"

# Build conda flags (omitted in container mode — mutually exclusive)
conda_flags=()
if [ -z "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    conda_flags+=(--use-conda --conda-frontend conda)
    [ -n "$SNAKEMAKE_CONDA_PREFIX" ] && conda_flags+=(--conda-prefix "$SNAKEMAKE_CONDA_PREFIX")
fi

# Build container flags
container_flags=()
if [ -n "$SNAKEMAKE_CONTAINER_FLAG" ]; then
    container_flags+=("$SNAKEMAKE_CONTAINER_FLAG")
    [ -n "$SNAKEMAKE_SIF_PREFIX" ]  && container_flags+=("$SNAKEMAKE_PREFIX_FLAG" "$SNAKEMAKE_SIF_PREFIX")
    [ -n "$SNAKEMAKE_BIND_PATHS" ]  && container_flags+=("$SNAKEMAKE_ARGS_FLAG" "$SNAKEMAKE_BIND_PATHS")
fi

snakemake --dryrun --cores 32 $snakemake_flags -s $snakefile \
    --printshellcmds --reason "${conda_flags[@]}" "${container_flags[@]}" $TARGETS
