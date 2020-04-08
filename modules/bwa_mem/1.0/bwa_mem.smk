#!/usr/bin/env snakemake
# Author: Helena Winata
# email/github: hwinata@bccrc.ca / whelena

##### SETUP #####
import os
import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "bwa_mem", 
    version = "1.0", 
    subdirs = ["inputs", "outputs"],
    req_references = ["bwa_index"]
)

localrules: 
    _bwa_mem_input,  
    _bwa_mem_all

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
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    log:
        bwa = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.bwa.stderr.log",
        samtools = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.samtools.stderr.log"
    params:
        bwa   = CFG["options"]["bwa_mem"],
        samtools  = CFG["options"]["samtools"],
        fasta  = lambda wildcards: config["reference"][wildcards.genome_build]["bwa_index"]
    conda:
        CFG["conda_envs"].get("bwa_mem") or "envs/bwa_mem.yaml"
    threads:
        CFG["threads"].get("bwa_mem") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("bwa_mem") or 10000
    shell:
        md.as_one_line("""
        bwa mem -t {threads} {params.bwa}
        {params.fasta} {input.fq} 
        2> {log.bwa} | samtools view {params.samtools}
        - > {output.bam} 2> {log.samtools}
        """)


rule _bwa_mem_all:
    input:
        fastqs = expand(rules._bwa_mem_run.output.bam, zip,
                        sample_id = list(CFG["samples"]['sample_id']), 
                        seq_type = list(CFG["samples"]['seq_type']),
                        genome_build = list(CFG["samples"]['genome_build']))

##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
