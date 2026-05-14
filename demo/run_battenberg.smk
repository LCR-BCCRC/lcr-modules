#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")
GENOME = op.filter_samples(SAMPLES, seq_type = "genome")


subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


configfile: "../modules/battenberg/1.2/config/default.yaml"
configfile: "genome_config.yaml"


config["lcr-modules"]["_shared"]["samples"] = GENOME


include: "../modules/battenberg/1.2/battenberg.smk"


rule all:
    input:
        rules._battenberg_all.input
