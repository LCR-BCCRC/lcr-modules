#!/usr/bin/env snakemake

##### PYTHON MODULES #####

import oncopipe as op
import pandas as pd
import os

##### SETUP VARIABLES #####

config["tool_names"] = "battenberg,cnaqc"
config["pipeline_name"] = "battenberg"  # adjust to match your battenberg runner's pipeline_name

include: "header_2.1.smk"

##### CONFIGURATION VALUES #####

configfile: "src/lcr-modules/modules/cnaqc/1.0/config/default.yaml"

config["lcr-modules"]["cnaqc"]["inputs"]["subclones"] = (
    "results/" + config["unix_group"] + "/battenberg-1.1"
    "/99-outputs/txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_subclones.txt"
)
config["lcr-modules"]["cnaqc"]["inputs"]["cellularity_ploidy"] = (
    "results/" + config["unix_group"] + "/battenberg-1.1"
    "/99-outputs/txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}_cellularity_ploidy.txt"
)
config["lcr-modules"]["cnaqc"]["inputs"]["sample_maf"] = (
    "results/" + config["unix_group"] + "/slms_3-1.0_vcf2maf-1.3"
    "/99-outputs/deblacklisted/maf"
    "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3.final.maf"
)

##### SAMPLE AND RUNS TABLES #####

# Use battenberg's runs table so CNAqc only runs on completed battenberg pairs
config["lcr-modules"]["_shared"]["samples"] = PIPELINE_SAMPLES["ALL"]["battenberg"]
config["lcr-modules"]["cnaqc"]["runs"] = PIPELINE_RUNS["ALL"]["battenberg"]

##### INCLUDE MODULES #####

include: "../../src/lcr-modules/modules/cnaqc/1.0/cnaqc.smk"

##### TARGETS #####

rule all:
    input:
        rules._cnaqc_all.input
