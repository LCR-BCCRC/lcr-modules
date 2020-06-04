#!/bin/bash
snakemake --keep-going --latency-wait 120 --jobs 100 --cluster-sync "srun -J {rule} --cpus-per-task {threads} --mem {resources.mem_mb}" \
--use-conda --conda-prefix /projects/clc/usr/anaconda/workflow_prototype/lcr-modules-envs/ all \
--rerun-incomplete
