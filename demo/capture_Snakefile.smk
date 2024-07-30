#!/usr/bin/env snakemake

'''
This Snakefile is made to run all the modules compatible with Capture workflow.
Compatibility of a workflow can be checked by referring to the pairing_config parameter present in a default.yaml file of that module.
'''
##### SETUP #####

import oncopipe as op

# filter sample table to use only capture seq_type
SAMPLES = op.load_samples("data/samples.tsv")
CAPTURE = op.filter_samples(SAMPLES, seq_type = "capture")


##### REFERENCE_FILES WORKFLOW #####


subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.5/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.5/config/default.yaml"


##### CONFIGURATION FILES #####


# Load module-specific configuration
configfile: "../modules/bwa_mem/1.1/config/default.yaml" # align samples
configfile: "../modules/utils/2.1/config/default.yaml"
configfile: "../modules/qc/1.0/config/default.yaml" # perform qc
configfile: "../modules/sequenza/1.4/config/default.yaml" # CNV calling
configfile: "../modules/battenberg/1.2/config/default.yaml"
configfile: "../modules/slms_3/1.0/config/default.yaml" # variant calling

# Load project-specific config, which includes the shared
# configuration and some module-specific config updates
configfile: "capture_config.yaml"


##### CONFIGURATION UPDATES #####


# Use all samples as a default sample list for each module
config["lcr-modules"]["_shared"]["samples"] = CAPTURE

##### MODULE SNAKEFILES #####


# Load module-specific snakefiles
include: "../modules/bwa_mem/1.1/bwa_mem.smk"
include: "../modules/utils/2.1/utils.smk"
include: "../modules/qc/1.0/qc.smk"
include: "../modules/sequenza/1.4/sequenza.smk"
include: "../modules/battenberg/1.2/battenberg.smk"
include: "../modules/slms_3/1.0/slms_3.smk"


##### TARGETS ######

rule all:
    input:
        rules._bwa_mem_all.input,
        rules._qc_all.input,
        rules._sequenza_all.input,
        rules._battenberg_all.input,
        rules._slms_3_all.input
        
