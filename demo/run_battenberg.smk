#!/usr/bin/env snakemake

'''
Example Snakefile for running the Battenberg module in isolation.

Must be run from the demo/ directory of lcr-modules. The symlink in
modules/battenberg/1.2/ is for reference only — paths will not work
if launched from the module directory.

    cd /path/to/lcr-modules/demo

    Dry run using conda:
    ./dry-run.sh run_battenberg.smk all "" runtime_config.conda.yaml

    Dry run using apptainer:
    ./dry-run.sh run_battenberg.smk all "" runtime_config.apptainer.yaml

    Local run using apptainer:
    ./run.sh run_battenberg.smk all "" runtime_config.apptainer.yaml

    SLURM cluster run using apptainer:
    ./snakemake.slurm.sh run_battenberg.smk all "" runtime_config.apptainer.yaml

    SLURM cluster run with additional snakemake arguments:
    ./snakemake.slurm.sh run_battenberg.smk all "--rerun-incomplete --keep-going" runtime_config.apptainer.yaml
'''

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")
CAPTURE = op.filter_samples(SAMPLES, seq_type = "capture")


subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


configfile: "../modules/battenberg/1.2/config/default.yaml"
configfile: "capture_config.yaml"


config["lcr-modules"]["_shared"]["samples"] = CAPTURE


include: "../modules/battenberg/1.2/battenberg.smk"


rule all:
    input:
        rules._battenberg_all.input
