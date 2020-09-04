#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Helena Winata, Bruno Grande
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op


##### RULES #####


# _utils_bam_sort: Sort a BAM file using coordinates
rule:
    input:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam"
    output:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.sort.bam"
    log:
        stdout = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_sort.stdout.log",
        stderr = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_sort.stderr.log"
    params:
        opts = CFG["options"].get("utils_bam_sort", ""),
        prefix = CFG["dirs"]["_parent"] + "{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CFG["threads"].get("utils_bam_sort", 12)
    resources: 
        mem_mb = CFG["mem_mb"].get("utils_bam_sort", 12000)
    shell:
        op.as_one_line("""
        samtools sort {params.opts} -@ {threads} -m $(({resources.mem_mb} / {threads} / 2))M
        -T {params.prefix} -o {output.bam} {input.bam} > {log.stdout} 2> {log.stderr}
        """)


# _utils_bam_markdups: Mark duplicates in a BAM file using Picard criteria
rule:
    input:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam"
    output:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.mdups.bam"
    log:
        stdout = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_markdups.stdout.log",
        stderr = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_markdups.stderr.log"
    params:
        opts = CFG["options"].get("utils_bam_markdups", ""),
        prefix = CFG["dirs"]["_parent"] + "{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("sambamba", "envs/sambamba-0.7.1.yaml")
    threads:
        CFG["threads"].get("utils_bam_markdups", 12)
    resources: 
        mem_mb = CFG["mem_mb"].get("utils_bam_markdups", 8000)
    shell:
        op.as_one_line("""
        sambamba markdup {params.opts} --nthreads {threads} 
        {input.bam} {output.bam} > {log.stdout} 2> {log.stderr}
        """)


# _utils_bam_index: Index a BAM file
rule:
    input:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam"
    output:
        bam = CFG["dirs"]["_parent"] + "{prefix}/{suffix}.bam.bai"
    log:
        stdout = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_index.stdout.log",
        stderr = CFG["logs"]["_parent"] + "{prefix}/{suffix}/bam_index.stderr.log"
    params:
        opts = CFG["options"].get("utils_bam_index", "-b"),
        prefix = CFG["dirs"]["_parent"] + "{prefix}/{suffix}"
    conda:
        CFG["conda_envs"].get("samtools", "envs/samtools-1.9.yaml")
    threads:
        CFG["threads"].get("utils_bam_index", 6)
    resources: 
        mem_mb = CFG["mem_mb"].get("utils_bam_index", 4000)
    shell:
        op.as_one_line("""
        samtools index {params.opts} -@ {threads} 
        {input.bam} > {log.stdout} 2> {log.stderr}
        """)
