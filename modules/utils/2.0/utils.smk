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

CONFIG = config["lcr-modules"]["utils"]
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"


def _get_log_dirs(wildcards):
    path = (wildcards.prefix).split(os.sep)
    root_dir = config["lcr-modules"]["_shared"]["root_output_dir"]

    # split root directory and remove 
    lroot = list(filter(None, root_dir.split(os.sep)))

    # get logs/launched...
    LOG = "logs/" + op._session.launched_fmt

    path.insert((len(lroot) + 1), LOG)
    path.append(wildcards.suffix)
    log_path = "/".join(path)
    print(log_path)

    logs = [log_path + "_bam_sort.stdout.log", log_path + "_bam_sort.stderr.log"]

    return logs


##### RULES #####
# _utils_bam_sort: Sort a BAM file using coordinates
rule:
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
        #logs = _get_log_dirs,
        opts = CONFIG["options"]["bam_sort"],
        prefix ="{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"]["samtools"]
    threads:
        CONFIG["threads"]["bam_sort"]
    resources: 
        mem_mb = CONFIG["mem_mb"]["bam_sort"]
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
rule:
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
        #logs = _get_log_dirs,
        opts = CONFIG["options"]["bam_markdups"],
        prefix = "{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"]["sambamba"]
    threads:
        CONFIG["threads"]["bam_markdups"]
    resources: 
        mem_mb = CONFIG["mem_mb"]["bam_markdups"]
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
rule:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.bam.bai"
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_index.stderr.log"
    params:
        #logs = _get_log_dirs,
        opts = CONFIG["options"]["bam_index"],
        prefix = "{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"]["samtools"]
    threads:
        CONFIG["threads"]["bam_index"]
    resources: 
        mem_mb = CONFIG["mem_mb"]["bam_index"]
    shell:
        op.as_one_line("""
        samtools index 
        {params.opts} 
        -@ {threads} 
        {input.bam} 
        > {log.stdout}
        2> {log.stderr}
        """)


rule: # create_interval_list from bed; default exomes
    input:
        bed = CONFIG["inputs"]["bed"],
        seq_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output: 
        intervals = "reference/exomes/{genome_build}/interval/{id}_intervals.txt"
    log: 
        "exomes/{genome_build}/interval/{id}_intervals.log"
    conda: 
        CONFIG["conda_envs"]["picard"]
    threads:
        CONFIG["threads"]["interval"]
    resources: 
        mem_mb = CONFIG["mem_mb"]["interval"]
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


del CONFIG

