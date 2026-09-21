#!/usr/bin/env snakemake

'''
Demo Snakefile: trust4/1.0

Runs TRUST4 BCR/TCR immune-repertoire reconstruction on bulk RNA-seq, using
modules/star/1.4 to align the demo's own mrna FASTQs first (TRUST4's own default
inputs.sample_bam points at star's final output). No licensing gate -- TRUST4 is
bioconda-installable, unlike mhc_hammer's Novoalign/HLA-HD/VEP.

Designed to run from the demo/ directory against the TCRBOA7 test dataset.
Usage (conda):      ./run.sh run_trust4.smk all "" runtime_config.conda.yaml
Usage (apptainer):  ./run.sh run_trust4.smk all "" runtime_config.apptainer.yaml
'''

##### PYTHON MODULES #####

import oncopipe as op


##### SAMPLES #####

SAMPLES = op.load_samples("data/samples.tsv")
MRNA    = op.filter_samples(SAMPLES, seq_type="mrna")


##### REFERENCE FILES WORKFLOW #####

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION #####

configfile: "../modules/utils/2.1/config/default.yaml"
configfile: "../modules/star/1.4/config/default.yaml"
configfile: "../modules/trust4/1.0/config/default.yaml"

configfile: "trust4_config.yaml"

# Use mrna samples only
config["lcr-modules"]["_shared"]["samples"] = MRNA


##### MODULES #####

include: "../modules/utils/2.1/utils.smk"
include: "../modules/star/1.4/star.smk"
include: "../modules/trust4/1.0/trust4.smk"


##### TARGETS #####

rule all:
    input:
        rules._trust4_all.input
