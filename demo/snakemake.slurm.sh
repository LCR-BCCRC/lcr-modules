#!/bin/bash

# Launches a gambl snakefile of your choice on numbers and sets the conda prefix

# Usage: ./snakemake.number.sh <snakefile.smk> <target_rule> <unix_group> <snakemake_flags>

# Example: ./snakemake.numbers.sh src/snakemake/manta.smk all gambl --rerun-incomplete

# snakefile.smk The snakefile you want to run
# target_rule: The name of one of the target rules specified in one of the included Snakefiles

# Set the ulimit to work with many concurrent jobs and files


snakefile=$1
target_rule=${2:-all}
snakemake_flags=$3

echo snakemake --use-conda --conda-prefix "/home/rmorin/miniconda3/envs" --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  --cluster-sync "srun -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 08:00:00" \
  -s $snakefile  $target_rule

sleep 5

snakemake --use-conda --conda-prefix "/home/rmorin/miniconda3/envs" --jobs 5000 \
  --latency-wait 120 $snakemake_flags --keep-going --default-resources mem_mb=2000 \
  --cluster-sync "srun -n 1 -N 1 -J {rule} --mem {resources.mem_mb} --cpus-per-task {threads} --time 08:00:00" \
  -s $snakefile  $target_rule 
