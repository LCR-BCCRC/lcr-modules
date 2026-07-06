#!/usr/bin/env snakemake

'''
Demo Snakefile: sage/1.2

Runs SAGE 5.1 somatic variant calling on paired tumour-normal WGS samples.
Uses the patched JAR bundled with the module (fixes for NPE on out-of-bounds
variants and fragment-sync performance). Container mode requires only a JRE
(eclipse-temurin:17-jre); the JAR is supplied from the module directory.

Designed to run from the demo/ directory against the TCRBOA7 test dataset.
Usage (conda):      ./run.sh run_sage.smk all "" runtime_config.conda.yaml
Usage (apptainer):  ./run.sh run_sage.smk all "" runtime_config.apptainer.yaml
'''

##### PYTHON MODULES #####

import oncopipe as op


##### SAMPLES #####

SAMPLES = op.load_samples("data/samples.tsv")
GENOME  = op.filter_samples(SAMPLES, seq_type="genome")


##### REFERENCE FILES WORKFLOW #####

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION #####

configfile: "../modules/sage/1.2/config/default.yaml"

# Point to the demo BAM files
config["lcr-modules"]["sage"]["inputs"]["sample_bam"] = "data/{sample_id}.bam"

# Use genome samples
config["lcr-modules"]["_shared"]["samples"] = GENOME


##### MODULE #####

include: "../modules/sage/1.2/sage.smk"


##### TARGETS #####

rule all:
    input:
        rules._sage_all.input
