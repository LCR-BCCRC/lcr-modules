#!/usr/bin/env snakemake

'''
It will only run through the workflow(eg. genome) sub directory.
'''
##### SETUP #####

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")
GENOME = op.filter_samples(SAMPLES, seq_type = "genome")


##### REFERENCE_FILES WORKFLOW #####


subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION FILES #####



# Load module-specific configuration

configfile: "../modules/slms_3/1.0/config/default.yaml"
configfile: "../modules/sage/1.0/config/default.yaml"
configfile: "../modules/manta/2.3/config/default.yaml"
#configfile: "../modules/gridss/1.1/config/default.yaml"
configfile: "../modules/strelka/1.1/config/default.yaml"
configfile: "../modules/lofreq/1.0/config/default.yaml"
configfile: "../modules/starfish/2.0/config/default.yaml"
configfile: "../modules/controlfreec/1.1/config/default.yaml"
configfile: "../modules/utils/2.1/config/default.yaml"
configfile: "../modules/picard_qc/1.0/config/default.yaml"
configfile: "../modules/bam2fastq/1.2/config/default.yaml"
#configfile: "../modules/vcf2maf/1.2/config/default.yaml"
configfile: "../modules/sequenza/1.4/config/default.yaml"
configfile: "../modules/bwa_mem/1.1/config/default.yaml"
configfile: "../modules/mixcr/1.1/config/default.yaml"
configfile: "../modules/battenberg/1.1/config/default.yaml"
configfile: "../modules/mutect2/2.0/config/default.yaml"
configfile: "../modules/varscan/1.1/config/default.yaml"
configfile: "../modules/liftover/1.2/config/default.yaml"
#configfile: "../modules/pathseq/1.0/config/default.yaml"
#configfile: "../modules/hmftools/1.0/config/default.yaml"

configfile: "../modules/hmftools/1.0/config/default.yaml"

# Load project-specific config, which includes the shared 
# configuration and some module-specific config updates
configfile: "genome_config.yaml"


##### CONFIGURATION UPDATES #####


# Use all samples as a default sample list for each module
config["lcr-modules"]["_shared"]["samples"] = GENOME

##### MODULE SNAKEFILES #####


# Load module-specific snakefiles
include: "../modules/slms_3/1.0/slms_3.smk"
include: "../modules/sage/1.0/sage.smk"
include: "../modules/manta/2.3/manta.smk"
include: "../modules/strelka/1.1/strelka.smk"
#include: "../modules/gridss/1.1/gridss.smk"
include: "../modules/lofreq/1.0/lofreq.smk"
include: "../modules/starfish/2.0/starfish.smk"
include: "../modules/controlfreec/1.1/controlfreec.smk"
include: "../modules/utils/2.1/utils.smk"
include: "../modules/picard_qc/1.0/picard_qc.smk"
#include: "../modules/vcf2maf/1.2/vcf2maf.smk"
include: "../modules/sequenza/1.4/sequenza.smk"
include: "../modules/bwa_mem/1.1/bwa_mem.smk"
include: "../modules/bam2fastq/1.2/bam2fastq.smk"
include: "../modules/sage/1.0/sage.smk"
include: "../modules/mixcr/1.1/mixcr.smk"
include: "../modules/battenberg/1.1/battenberg.smk"
include: "../modules/mutect2/2.0/mutect2.smk"
include: "../modules/varscan/1.1/varscan.smk"
include: "../modules/liftover/1.2/liftover.smk"
#include: "../modules/pathseq/1.0/pathseq.smk"
#include: "../modules/hmftools/1.0/hmftools.smk"



##### TARGETS ######

rule all:
    input:
        rules._sage_all.input,
        rules._picard_qc_all.input,
        rules._bam2fastq_all.input,
        rules._sequenza_all.input,
        rules._bwa_mem_all.input,
        rules._controlfreec_all.input,
        #rules._vcf2maf_all.input
        rules._slms_3_all.input,
        rules._manta_all.input,  
        rules._lofreq_all.input,
        rules._strelka_all.input,
        #rules._gridss_all.input,
        rules._starfish_all.input,
        rules._mixcr_all.input,
        rules._mutect2_all.input,
        rules._varscan_all.input,
        rules._liftover_all.input,
        rules._battenberg_all.input,
        #rules._pathseq_all.input,
        #rules._hmftools_all.input

