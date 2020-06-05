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

CFG = config["lcr-modules"]["utils"]

def _get_log_dirs(wildcards):
    path = (wildcards.prefix).split(os.sep)
    root_dir = config["lcr-modules"]["_shared"]["root_output_dir"]

    # split root directory and remove 
    lroot = list(filter(None, root_dir.split(os.sep)))

    # get logs/launched...
    LOG = "logs/" + op._session.launched_fmt

    path.insert((len(lroot) + 1), LOG)
    log_path = "/".join(path)
    print(log_path)

    logs = [log_path + "bam_sort.stdout.log", log_path + "bam_sort.stderr.log"]

    return logs


##### RULES #####
# _utils_bam_sort: Sort a BAM file using coordinates
rule:
    input:
        bam = "{prefix}/{suffix}.bam"
    output:
        bam = "{prefix}/{suffix}.sort.bam"
    params:
        logs = _get_log_dirs,
        opts = CFG["options"].get("bam_sort", ""),
        prefix ="{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CFG["threads"].get("bam_sort", 12)
    resources: 
        mem_mb = CFG["mem_mb"].get("bam_sort", 12000)
    shell:
        op.as_one_line("""
        samtools sort {params.opts} -@ {threads} -m $(({resources.mem_mb} / {threads}))M
        -T {params.prefix} -o {output.bam} {input.bam} 
        > {params.logs}[0]
        2> {params.logs}[1]
        """)


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule:
    input:
        bam = "{prefix}/{suffix}.bam"
    output:
        bam = "{prefix}/{suffix}.mdups.bam"
    params:
        logs = _get_log_dirs,
        opts = CFG["options"].get("bam_markdups", ""),
        prefix = "{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("sambamba", "envs/sambamba-0.7.1.yaml")
    threads:
        CFG["threads"].get("bam_markdups", 12)
    resources: 
        mem_mb = CFG["mem_mb"].get("bam_markdups", 8000)
    shell:
        op.as_one_line("""
        sambamba markdup {params.opts} --nthreads {threads} 
        {input.bam} {output.bam} 
        > {params.logs}[0]
        2> {params.logs}[1]
        """)


# _utils_bam_index: Index a BAM file
rule:
    input:
        bam = "{prefix}/{suffix}.bam"
    output:
        bam = "{prefix}/{suffix}.bam.bai"
    params:
        logs = _get_log_dirs,
        opts = CFG["options"].get("bam_index", "-b"),
        prefix = "{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CFG["threads"].get("bam_index", 6)
    resources: 
        mem_mb = CFG["mem_mb"].get("bam_index", 4000)
    shell:
        op.as_one_line("""
        samtools index {params.opts} -@ {threads} 
        {input.bam} 
        > {params.logs}[0]
        2> {params.logs}[1]
        """)

del CFG