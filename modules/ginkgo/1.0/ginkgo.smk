#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "ginkgo", 
    version = "1.0",
    # TODO: If applicable, add other output subdirectories
    subdirs = ["inputs", "bed", "ginkgo", "outputs"],
    # TODO: Replace "genome_fasta" with actual reference requirements
    req_references = ["genome_fasta"] 
)

# Define rules to be run locally when using a compute cluster.
# TODO: Replace with actual rule names
localrules: 
    _ginkgo_input_bed,
    _ginkgo_output,
    _ginkgo_all


##### RULES #####
# in single cell modules, sample_id and lib_id refers to patient_id and sample_id respectively on the datasheet
# for ginkgo purposes, each library = a cell

def get_input_bam(wildcards):


rule _ginkgo_input_bam:
    input:
        bam = CFG["inputs"]["sample_bed"]
        "results/bwa_mem-1.0/99-outputs/scWGS--{genome_build}/{lib_id}.bam"
    output:
        bam = CFG["dirs"]["inputs"] + "{genome_build}/{sample_id}/{lib_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _ginkgo_bam2bed:
    input:
        bam = CFG["dirs"]["inputs"] + "{genome_build}/{sample_id}/{lib_id}.bam"
    output:
        bed = CFG["dirs"]["bed"] + "/{sample}/{lib}.bed.gz"
    conda:
        "/projects/clc/usr/anaconda/workflow_prototype/envs/484bc115.yaml"
    log:
        "logs/{sample}/{lib}_bedtools.log"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        "bedtools bamtobed -i {input} > >(gzip > {output}) 2> {log}"


rule _ginkgo_link_bed_to_bins:
    input:
        bed = rules._ginkgo_bam2bed
    output:
        bed = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/{lib_id}.bed.gz"
    run:
        md.symlink(input.bed, output.bed)


def get_bed_files(wildcards):
    sub_df = md.filter_samples(CFG["samples", patient_id = wildcards.sample_id])
    libs = list(sub_df["sample_id"])
    bed = expand("{dir}{{genome_build}}_bin{{bin}}/{sample_id}/{lib_id}.bed.gz", dir = CFG["dirs"]["ginkgo"], sample_id = wildcards.sample_id, lib_id = libs)
    return bed


rule _ginkgo_run:
    input:
        bed = get_bed_files
    output:
        config["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.done"
    log: 
        stdout = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stderr.log"
    params:
        ginkgo = CFG["software"],
        bedDir = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}",
        genome = {genome_build}
        binning = CFG["options"]["binMeth"],
        clustDist = CFG["options"]["distMeth"],
        clustLink = CFG["options"]["clustMeth"],
        opts = CFG["options"]["flags"]
    threads: 2
    resources:
        mem_mb = 8000
    shell:
        as_one_line("""
        {params.ginkgo} 
        --input {params.bedDir} 
        --genome {params.genome} 
        --binning {params.binning}
        {params.clustDist} {params.clustLink}
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        && touch {output}
        """)


rule _ginkgo_output:
    input:
        sc = config["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/SegCopy"
    output:
        sc = CFG["dirs"]["outputs"] + "{genome_build}_bin{bin}/{sample_id}_SegCopy"
    run:
        md.symlink(input.sc, output.sc)

sub_df = CFG["samples"][["patient_id", "genome_build"]].drop_duplicates()

rule _ginkgo_all:
    input:
        expand(expand("{{dir}}{genome_build}_bin{{bin}}/{sample_id}_SeqCopy", zip,
               genome_build = sub_df["genome_build"],
               sample_id = sub_df["patient_id"]),
            dir = CFG["dirs"]["outputs"]
            bin = CFG["inputs"]["bins"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
