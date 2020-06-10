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

print("utilities module is included")
CONFIG = config["lcr-modules"]["utils"]
PATH = ".*\/(" + "|".join(CONFIG["paired_modules"]) + ").*\/"

'''

for m in CONFIG["paired_modules"]:
    print (config["lcr-modules"])
    pdir = config["lcr-modules"][m]["dirs"]["_parent"]
    print(p_dir)
'''

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

wildcard_constraints:
    prefix = "[0-9]{2}-.*"

##### RULES #####
# _utils_bam_sort: Sort a BAM file using coordinates
rule:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.sort.bam"
    priority: 5
    params:
        logs = _get_log_dirs,
        opts = CONFIG["options"].get("bam_sort", ""),
        prefix ="{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CONFIG["threads"].get("bam_sort", 12)
    resources: 
        mem_mb = CONFIG["mem_mb"].get("bam_sort", 12000)
    shell:
        op.as_one_line("""
        samtools sort {params.opts}
        -@ {threads} -m $(({resources.mem_mb} / {threads}))M
        -T {params.prefix} -o {output.bam} {input.bam} 
        > {params.logs[0]}
        2> {params.logs[1]}
        """)


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.mdups.bam"
    priority: 5
    params:
        logs = _get_log_dirs,
        opts = CONFIG["options"].get("bam_markdups", ""),
        prefix = "{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"].get("sambamba", "envs/sambamba-0.7.1.yaml")
    threads:
        CONFIG["threads"].get("bam_markdups", 12)
    resources: 
        mem_mb = CONFIG["mem_mb"].get("bam_markdups", 8000)
    shell:
        op.as_one_line("""
        sambamba markdup {params.opts} 
        --nthreads {threads} 
        {input.bam} {output.bam} 
        > {params.logs[0]}
        2> {params.logs[1]}
        """)


# _utils_bam_index: Index a BAM file
rule:
    input:
        bam = "{out_dir}/{prefix}/{suffix}.bam"
    output:
        bam = "{out_dir}/{prefix}/{suffix}.bam.bai"
    priority: 5
    params:
        logs = _get_log_dirs,
        opts = CONFIG["options"].get("bam_index", "-b"),
        prefix = "{out_dir}/{prefix}/{suffix}"
    conda:
        CONFIG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CONFIG["threads"].get("bam_index", 6)
    resources: 
        mem_mb = CONFIG["mem_mb"].get("bam_index", 4000)
    shell:
        op.as_one_line("""
        samtools index {params.opts} -@ {threads} 
        {input.bam} 
        > {params.logs[0]}
        2> {params.logs[1]}
        """)

del CONFIG
