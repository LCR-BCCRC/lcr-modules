#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Helena Winata, Bruno Grande
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####

import sys, os
from os.path import join

# Import package with useful functions for developing analysis modules
import oncopipe as op

_UTILS = config["lcr-modules"]["utils"]
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"


##### RULES #####
# _utils_bam_sort: Sort a BAM file using coordinates
rule _utils_bam_sort:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.sort.bam"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_sort.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_sort.stderr.log"
    wildcard_constraints:
        prefix = ".*(sort).*"
    params:
        opts = _UTILS["options"]["bam_sort"],
        prefix ="{out_dir}/{prefix}/{suffix}"
    conda:
        _UTILS["conda_envs"]["samtools"]
    threads:
        _UTILS["threads"]["bam_sort"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_sort"]
    shell:
        op.as_one_line("""
        samtools sort 
        -@ {threads} -m $(({resources.mem_mb} / {threads}))M
        {params.opts}
        -T {params.prefix} -o {output.bam} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule _utils_bam_markdups:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.mdups.bam"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_mark_dups.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_mark_dups.stderr.log"
    wildcard_constraints:
        prefix = ".*(mark_dups).*"
    params:
        opts = _UTILS["options"]["bam_markdups"]
    conda:
        _UTILS["conda_envs"]["sambamba"]
    threads:
        _UTILS["threads"]["bam_markdups"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_markdups"]
    shell:
        op.as_one_line("""
        sambamba markdup 
        {params.opts} 
        --nthreads {threads} 
        {input.bam} {output.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


# _utils_bam_index: Index a BAM file
rule _utils_bam_index:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.bam.bai"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stderr.log"
    params:
<<<<<<< HEAD
        opts = _UTILS["options"]["bam_index"],
=======
        opts = _UTILS["options"]["bam_index"]
>>>>>>> 664903ee9d6f0c4c9e3ba2f20dc7edee662b572e
    conda:
        _UTILS["conda_envs"]["samtools"]
    threads:
        _UTILS["threads"]["bam_index"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["bam_index"]
    shell:
        op.as_one_line("""
        samtools index 
        {params.opts} 
        -@ {threads} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)
<<<<<<< HEAD


rule _utils_create_intervals: # create_interval_list from bed; default exomes
    input:
        bed = lambda w: _UTILS["inputs"]["bed"][w.genome_build],
        seq_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        intervals = "reference/exomes/{genome_build}/interval/{id}_intervals.txt"
    log: 
        "exomes/{genome_build}/interval/{id}_intervals.log"
    conda: 
        _UTILS["conda_envs"]["picard"]
    threads:
        _UTILS["threads"]["interval"]
    resources: 
        mem_mb = _UTILS["mem_mb"]["interval"]
    shell:
        op.as_one_line("""
        picard BedToIntervalList
        I={input.bed}
        O={output.intervals}
        SD={input.seq_dict}
        2> {log}
        &&
        chmod a-w {output.intervals}
        """)



=======
>>>>>>> 664903ee9d6f0c4c9e3ba2f20dc7edee662b572e
