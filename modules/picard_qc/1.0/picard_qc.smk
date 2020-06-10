#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["picard_qc"]`
CFG = op.setup_module(
    name = "picard_qc",
    version = "1.0",
    subdirectories = ["inputs", "metrics", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _picard_qc_input_bam,
    _picard_qc_rrna_int,
    _picard_qc_output_txt,
    _picard_qc_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _picard_qc_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _picard_qc_alignment_summary:
    input:
        bam = rules._picard_qc_input_bam.output.bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        summary = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary_metrics",
        insert_size = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/insert_size_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/alignment_summary.stderr.log"
    params:
        opts = CFG["options"]["alignment_summary"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["alignment_summary"]
    resources:
        mem_mb = CFG["mem_mb"]["alignment_summary"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectMultipleMetrics {params.opts} 
        I={input.bam} O=$( dirname {output.summary}) R={input.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_hs_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        intervals = reference_files(CFG["inputs"]["intervals"])
    output:
        hs = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics",
        intervals = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/interval_hs_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/hs_metrics.stderr.log"
    params:
        opts = CFG["options"]["hs_metrics"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["hs_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["hs_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectHsMetrics {params.opts} 
        I={input.bam} O={output.hs} PER_TARGET_COVERAGE={output.intervals} 
        R={input.fasta} TI={input.intervals} BI={input.intervals} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_rrna_int:
    input:
        bam = rules._picard_qc_input_bam.output.bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        rrna_int = CFG["inputs"]["rrna_intervals"]
    output:
        rrna_int = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rrna_int_list"
    log:
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rrna_int.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        op.as_one_line("""
        samtools view -H {input.bam} | 
        grep -v "@RG" | cat - {input.rrna_int} 
        > {output.rrna_int}
        2> {log.stderr}
        """)


rule _picard_qc_rnaseq_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.bam,
        rrna_int = rules._picard_qc_rrna_int.output.rrna_int,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        refFlat = CFG["inputs"]["refFlat"]
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/rnaseq_metrics.stderr.log"
    params:
        opts = CFG["options"]["rnaseq_metrics"]["base"],
        strand = op.switch_on_column("strand", CFG["samples"], CFG["options"]["rnaseq_metrics"]["strand"], match_on = "sample")
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["rnaseq_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["rnaseq_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectRnaSeqMetrics {params.opts} 
        I={input.bam} O={output.metrics} 
        R={input.fasta} 
        STRAND_SPECIFICITY={params.strand} 
        REF_FLAT={input.refFlat}
        RIBOSOMAL_INTERVALS={input.rrna_int}
        CHART_OUTPUT={output.metrics}.pdf 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_wgs_metrics:
    input:
        bam = rules._picard_qc_input_bam.output.bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics"
    log:
        stdout = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics.stdout.log",
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/wgs_metrics.stderr.log"
    params:
        opts = CFG["options"]["wgs_metrics"],
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["wgs_metrics"]
    resources:
        mem_mb = CFG["mem_mb"]["wgs_metrics"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectWgsMetrics {params.opts} 
        I={input.bam} O={output.metrics} 
        R={input.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_flagstats:
    input:
        bam = rules._picard_qc_input_bam.output.bam
    output:
        flagstats = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/flagstats"
    log:
        stderr = CFG["logs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/flagstats.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["flagstats"]
    resources:
        mem_mb = CFG["mem_mb"]["flagstats"]
    shell:
        op.as_one_line("""
        samtools flagstat -@ {threads}
        {input.bam} > {output.flagstats} 
        2> {log.stderr}
        """)


rule _picard_qc_output:
    input:
        metrics = CFG["dirs"]["metrics"] + "{seq_type}--{genome_build}/{sample_id}/{qc_type}"
    output:
        metrics = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{qc_type}/{sample_id}.{qc_type}"
    wildcard_constraints:
        qc_type = ".*_metrics|flagstats"
    run:
        op.relative_symlink(input.metrics, output.metrics)


def _get_picard_qc_files(wildcards):
    # base metrics to calculate for all seq_type
    base = ["alignment_summary_metrics", "insert_size_metrics", "flagstats"]
    # m_base = ["alignment_summary_metrics", "insert_size_metrics"]
    # additional seq-type-specific metrics
    if wildcards.seq_type == 'mrna':
        metrics = base + ["rnaseq_metrics"]
    #    m_metrics = mbase + ["rnaseq_metrics"]
    elif wildcards.seq_type == 'genome':
        metrics = base + ["wgs_metrics"]
    #    m_metrics = mbase + ["wgs_metrics"]
    elif wildcards.seq_type == 'capture':
        metrics = base + ["hs_metrics", "interval_hs_metrics"]
    #    m_metrics = mbase + ["hs_metrics"]
    else:
        metrics = base
    #    m_metrics = m_base

    targets = expand(rules._picard_qc_output.output.metrics, **wildcards, qc_type = metrics)

    # merged_targets = expand("{dir}{seq_type}--{genome_build}/all.{qc_type}.txt", dir = CFG["dirs"]["merged_metrics"], **wildcards, qc_type = m_metrics)

    return targets

 
rule _picard_qc_all_dispatch:
    input:
        _get_picard_qc_files
    output:
        touch(CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.dispatched")


rule _picard_qc_all:
    input:
        expand(rules._picard_qc_all_dispatch.output, zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
