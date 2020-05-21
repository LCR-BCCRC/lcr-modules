#!/usr/bin/env snakemake

##### SETUP #####

import os
import gzip
import modutils as md

CFG = md.setup_module(
    config = config, 
    name = "fastqc", 
    version = "1.0",
    subdirs = ["inputs", "fastqc", "outputs"],
    req_references = []
)

localrules: _fastqc_input, _fastqc_output, _fastqc_all

##### RULES #####

rule _fastqc_input:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _fastqc_run:
    input:
        bam = rules._fastqc_input.output.bam
    output:
        qc = CFG["dirs"]["fastqc"] + "{seq_type}--{genome_build}/{sample_id}_fastqc.zip"
    log:
        stdout = CFG["logs"]["fastqc"] + "{seq_type}--{genome_build}/{sample_id}/fastqc.stdout.log",
        stderr = CFG["logs"]["fastqc"] + "{seq_type}--{genome_build}/{sample_id}/fastqc.stderr.log"
    params:
        opts = CFG["options"]["fastqc"],
        dir = CFG["dirs"]["fastqc"] + "{seq_type}--{genome_build}"
    conda:
        CFG["conda_envs"].get("fastqc") or "envs/fastqc-0.11.9.yaml"
    threads:
        CFG["threads"].get("fastqc") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("fastqc") or 6000
    shell:
        md.as_one_line("""
        fastqc {params.opts} 
        -t {threads} 
        -o {params.dir} {input.bam} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _fastqc_output:
    input:
        qc = rules._fastqc_run.output.qc
    output:
        qc = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}_fastqc.zip"
    run:
        md.symlink(input.qc, output.qc)


rule _fastqc_all:
    input:
        qc = expand(rules._fastqc_output.output.qc, zip,
                    seq_type = list(CFG["samples"]["seq_type"]),
                    genome_build = list(CFG["samples"]["genome_build"]),
                    sample_id = list(CFG["samples"]["sample_id"]))

##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
