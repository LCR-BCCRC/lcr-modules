#!/usr/bin/env snakemake
# -----------------------------------------------------------------------------
# Author:        Helena Winata
# email/github:  hwinata@bccrc.ca / whelena
# -----------------------------------------------------------------------------
# Input:         {sample_id}.1.fastq         
#                {sample_id}.2.fastq         
#                                           
# Output:        {sample_id}.bam  
# 
# Required references:  genome.fa + indexes (genome.fa.fai, genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.pas, genome.fa.sa)
#
# Purpose: Align paired-end fastq reads (generated using Illumina) using a reference genome to create an aligned BAM file                     
# -----------------------------------------------------------------------------

##### SETUP #####
import os
import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "bwa_mem", 
    version = "1.0", 
    subdirs = ["inputs", "bam_util", "outputs"],
    req_references = ["bwa_index"],
    scratch_subdirs = ["bam_util"]
)

localrules: 
    _bwa_mem_input,
    _bwa_mem_output,
    _bwa_mem_all

# utility rules
include: "bam_util.smk"

##### RULES #####

rule _bwa_mem_input:
    input:
        fq = CFG["inputs"].get("fastq")
    output:
        fq = expand("{fqDIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{num}.fastq.gz", fqDIR = CFG["dirs"]["inputs"], num = ["1", "2"]) 
    run:
        md.symlink(input.fq[0], output.fq[0])
        md.symlink(input.fq[1], output.fq[1])


rule _bwa_mem_run:
    input:
        fq = rules._bwa_mem_input.output.fq
    output:
        sam = pipe(CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}.sam")
    log:
        bwa = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}.bwa.stderr.log",
    params:
        bwa   = CFG["options"]["bwa_mem"],
        fai  = md.get_reference(CFG, "bwa_index")
    conda:
        CFG["conda_envs"].get("bwa_mem") or "envs/bwa-0.7.17.yaml"
    threads:
        CFG["threads"].get("bwa_mem") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("bwa_mem") or 10000
    shell:
        md.as_one_line("""
        bwa mem -t {threads} {params.bwa}
        {params.fai} {input.fq} 
        > {output.sam} 2> {log.bwa}
        """)


rule _bwa_mem_samtools:
    input:
        sam = rules._bwa_mem_run.output.sam
    output:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    log:
        samtools = CFG["logs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}.samtools.stderr.log"
    params:
        samtools  = CFG["options"]["samtools"]
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("samtools") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("samtools") or 5000
    shell:
        md.as_one_line("""
        samtools view {params.samtools}
        {input.sam} > {output.bam} 2> {log.samtools}
        """)


rule _bwa_mem_output:
    input:
        bam = CFG["dirs"]["bam_util"] + "{seq_type}--{genome_build}/{sample_id}" + CFG["suffix"] + ".bam"
    output:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    shell:
        "ln {input} {output}"


rule _bwa_mem_all:
    input:
        bai = expand(rules._bam_util_index.output.bai, zip,
                    seq_type = list(CFG["samples"]['seq_type']),
                    genome_build = list(CFG["samples"]['genome_build']),
                    sample_id = list(CFG["samples"]['sample_id']))
'''
       fastqs = expand(rules._bwa_mem_run.output.bam, zip,
                        sample_id = list(CFG["samples"]['sample_id']), 
                        seq_type = list(CFG["samples"]['seq_type']),
                        genome_build = list(CFG["samples"]['genome_build']))
'''
##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
