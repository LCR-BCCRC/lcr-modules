#!/usr/bin/env snakemake

'''
Example Snakefile for running the Battenberg module in isolation.

This file must be run from the demo/ directory of lcr-modules:

    cd /path/to/lcr-modules/demo
    snakemake --snakefile run_battenberg.smk --use-singularity [other options]

The symlink in modules/battenberg/1.2/ is provided for reference only.
All paths in this file are relative to demo/ and will not work if the
file is launched from the module directory.
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
