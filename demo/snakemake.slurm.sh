#!/bin/bash

# Launches a snakefile of your choice on a Slurm cluster and sets the conda prefix

# Usage: ./snakemake.slurm.sh <snakefile.smk> <target_rule> <conda_prefix> <snakemake_flags>

# Example: ./snakemake.slurm.sh src/snakemake/manta.smk all "--rerun-incomplete" ~/minconda3/envs

# snakefile.smk The snakefile you want to run
# target_rule: The name of one of the target rules specified in one of the included Snakefiles
# snakemake_flags: Any additional Snakemake options (e.g. "-np" for a dry run)
# conda_prefix: The location where conda packages should be installed

snakefile=$1
target_rule=${2:-all}
snakemake_flags=$3
conda_prefix=${4:-$CONDA_PREFIX}


echo snakemake --use-conda --conda-prefix $conda_prefix --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  --cluster-sync "srun -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 24:00:00" \
  -s $snakefile  $target_rule

sleep 5

snakemake --use-conda --conda-prefix $conda_prefix --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  --cluster-sync "srun -p upgrade -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 08:00:00" \
  -s $snakefile  $target_rule 
