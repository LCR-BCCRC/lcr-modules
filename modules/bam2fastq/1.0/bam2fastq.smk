#!/usr/bin/env snakemake
# Author: Helena Winata
# email/github: hwinata@bccrc.ca / whelena


##### SETUP #####

import os
import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "bam2fastq", 
    version = "1.0",
    subdirs = ["inputs", "outputs"],
    req_references = ["genome_fasta"]
)

localrules: 
    _bam2fastq_input, 
    _bam2fastq_all

##### RULES #####

rule _bam2fastq_input:
    input:
        bam = CFG["inputs"].get("sample_bam")
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _bam2fastq_run:
    input:
        bam = rules._bam2fastq_input.output.bam
    output:
        fq = expand("{fqDIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{readNum}.fastq.gz", 
                     fqDIR = CFG["dirs"]["outputs"], readNum = ['1', '2'])
    log:
        stdout = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}/bam2fastq.stdout.log",
        stderr = CFG["logs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}/bam2fastq.stderr.log"
    params:
        opts   = CFG["options"]["picard"],
        fasta  = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("bam2fastq") or "envs/bam2fastq.yaml"
    threads:
        CFG["threads"].get("bam2fastq") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("bam2fastq") or 5000
    shell:
        md.as_one_line("""
        picard -Xmx{resources.mem_mb}m SamToFastq {params.opts}
        I={input.bam} FASTQ=>(gzip > {output.fq[0]}) SECOND_END_FASTQ=>(gzip > {output.fq[1]}) 
        > {log.stdout} &> {log.stderr}
        """)


rule _bam2fastq_all:
    input:
        fastqs = expand(rules._bam2fastq_run.output.fq, seq_type = CFG["runs"]["tumour_seq_type"], genome_build = CFG["runs"]["tumour_genome_build"], sample_id = list(CFG["samples"]['sample_id']))
        
        #expand("{dir}{seq_type}--{genome_build}/{sample_id}.{readNum}.fastq.gz", dir = CFG["dirs"]["outputs"], seq_type = CFG["runs"]["tumour_seq_type"], genome_build = CFG["runs"]["tumour_genome_build"], sample_id = list(CFG["samples"]['sample_id']), readNum = ['1', '2']) 


##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
