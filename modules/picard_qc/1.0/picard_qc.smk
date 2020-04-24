#!/usr/bin/env snakemake

##### SETUP #####
import os
import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "picard_qc", 
    version = "1.0",
    subdirs = ["inputs", "alignment_summary", "hs_metrics", "rnaseq_metrics", "wgs_metrics", "flagstats", "outputs"],
    req_references = ["genome_fasta"]
)

localrules: 
    _picard_qc_input, 
    _picard_qc_rrna_int,
    _picard_qc_output, 
    _picard_qc_all


##### RULES #####

rule _picard_qc_input:
    input:
        bam = CFG["inputs"].get("sample_bam")
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _picard_qc_alignment_summary:
    input:
        bam = rules._picard_qc_input.output.bam
    output:
        summary = CFG["dirs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.alignment_summary_metrics",
        insert_size = CFG["dirs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.insert_size_metrics"
    log:
        stdout = CFG["logs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.alignment_sum.stdout.log",
        stderr = CFG["logs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.alignment_sum.stderr.log"
    params:
        fasta  = md.get_reference(CFG, "genome_fasta"),
        opts = CFG["options"]["alignment_summary"],
        prefix = CFG["dirs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}"
    conda: 
        CFG["conda_envs"].get("picard") or "envs/picard-2.21.7.yaml"
    threads:
        CFG["threads"].get("alignment_summary") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("alignment_summary") or 6000
    shell:
        md.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectMultipleMetrics {params.opts} 
        I={input.bam} O={params.prefix} R={params.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_hs_metrics:
    input:
        bam = rules._picard_qc_input.output.bam
    output:
        hs = CFG["dirs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.hs_metrics",
        int = CFG["dirs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.interval_hs_metrics"
    log:
        stdout = CFG["logs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.hs_metrics.stdout.log",
        stderr = CFG["logs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.hs_metrics.stderr.log"
    params:
        opts  = CFG["options"]["hs_metrics"],
        fasta = md.get_reference(CFG, "genome_fasta"),
        int =  md.get_reference(CFG, "intervals")
    conda:
        CFG["conda_envs"].get("picard") or "envs/picard-2.21.7.yaml"
    threads:
        CFG["threads"].get("hs_metrics") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("hs_metrics") or 5000
    shell:
        md.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectHsMetrics {params.opts} 
        I={input.bam} O={output.hs}PER_TARGET_COVERAGE={output.int} 
        R={params.fasta} TI={params.int} BI={params.int} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _picard_qc_rrna_int:
    input:
        bam = rules._picard_qc_input.output.bam
    output:
        int = CFG["dirs"]["rnaseq_metrics"]+ "{seq_type}--{genome_build}/{sample_id}.rrna_int_list"
    log:
        CFG["logs"]["rnaseq_metrics"] + "{seq_type}--{genome_build}/{sample_id}.rrna_int.stderr.log"
    params: 
        rRNA = md.get_reference(CFG, "rRNA_int")
    conda: 
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    shell:
        'samtools view -H {input.bam} | grep -v "@RG" | cat - {params.rRNA}'


rule _picard_qc_rnaseq_metrics:
    input:
        bam = rules._picard_qc_input.output.bam,
        rRNA_int = rules._picard_qc_rrna_int.output.int
    output:
        metrics = CFG["dirs"]["rnaseq_metrics"]+ "{seq_type}--{genome_build}/{sample_id}.collect_rnaseq_metrics"
    log:
        stdout = CFG["logs"]["rnaseq_metrics"] + "{seq_type}--{genome_build}/{sample_id}.rnaseq_metrics.stdout.log",
        stderr = CFG["logs"]["rnaseq_metrics"] + "{seq_type}--{genome_build}/{sample_id}.rnaseq_metrics.stderr.log"
    params:
        opts  = CFG["options"]["rnaseq_metrics"]["base"],
        strand = md.switch_on_column("strand", CFG["samples"], CFG["options"]["rnaseq_metrics"]["strand"]),
        fasta = md.get_reference(CFG, "genome_fasta"),
        refFlat = md.get_reference(CFG, "refFlat")
    conda: 
        CFG["conda_envs"].get("picard") or "envs/picard-2.21.7.yaml"
    threads:
        CFG["threads"].get("rnaseq_metrics") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("rnaseq_metrics") or 5000
    shell:
        md.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectRnaSeqMetrics {params.opts} 
        I={input.bam} O={output.metrics} R={params.fasta} 
        STRAND_SPECIFICITY={params.strand} 
        REF_FLAT={params.refFlat}
        RIBOSOMAL_INTERVALS={input.rRNA_int}
        CHART_OUTPUT={output.metrics}.pdf 
        > {log.stdout} 2> {log.stderr}
        """)


rule _picard_qc_wgs_metrics:
    input:
        bam = rules._picard_qc_input.output.bam
    output:
        metrics = CFG["dirs"]["wgs_metrics"]+ "{seq_type}--{genome_build}/{sample_id}.collect_wgs_metrics"
    log:
        stdout = CFG["logs"]["wgs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.wgs_metrics.stdout.log",
        stderr = CFG["logs"]["wgs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.wgs_metrics.stderr.log"
    params:
        opts  = CFG["options"]["wgs_metrics"],
        fasta = md.get_reference(CFG, "genome_fasta")
    conda: 
        CFG["conda_envs"].get("picard") or "envs/picard-2.21.7.yaml"
    threads:
        CFG["threads"].get("wgs_metrics") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("wgs_metrics") or 8000
    shell:
        md.as_one_line("""
        picard -Xmx{resources.mem_mb}m CollectWgsMetrics {params.opts} 
        I={input.bam} O={output.metrics} R={params.fasta} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _picard_qc_flagstats:
    input:
        bam = rules._picard_qc_input.output.bam
    output:
        flagstats = CFG["dirs"]["flagstats"]+ "{seq_type}--{genome_build}/{sample_id}.flagstats"
    log:
        CFG["logs"]["flagstats"] + "{seq_type}--{genome_build}/{sample_id}.flagstats.stderr.log"
    conda: 
        CFG["conda_envs"].get("samtools") or "envs/samtools-1.9.yaml"
    threads:
        CFG["threads"].get("flagstats") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("flagstats") or 8000
    shell:
        "samtools flagstat {input.bam} > {output.flagstats} 2> {log}"


rule _picard_qc_output:
    input:
        metrics = md.switch_on_wildcard("qc_type", )"{seq_type}--{genome_build}/{sample_id}.{qc_type}"
    output:
        metrics = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.{qc_type}"
    run:
        md.symlink(input.metrics, output.metrics)

'''
outputs:
        summary = CFG["dirs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.alignment_summary_metrics",
        insert_size = CFG["dirs"]["alignment_summary"] + "{seq_type}--{genome_build}/{sample_id}.insert_size_metrics"
        flagstats = CFG["dirs"]["flagstats"]+ "{seq_type}--{genome_build}/{sample_id}.flagstats"

        hs = CFG["dirs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.hs_metrics",
        int = CFG["dirs"]["hs_metrics"] + "{seq_type}--{genome_build}/{sample_id}.interval_hs_metrics"
        metrics = CFG["dirs"]["rnaseq_metrics"]+ "{seq_type}--{genome_build}/{sample_id}.collect_rnaseq_metrics"
        metrics = CFG["dirs"]["wgs_metrics"]+ "{seq_type}--{genome_build}/{sample_id}.collect_wgs_metrics"

def _get_picard_qc_files(wildcards)
    DIR = CFG["dirs"]["outputs"] + "{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.sample_id}."
    base = ["alignment_summary_metrics", "insert_size_metrics", "flagstats"]
'''

rule _picard_qc_all_dispatch:
    input:
        base = expand("{DIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{qc_type}", DIR = CFG["dirs"]["outputs"], qc_type = ["alignment_summary_metrics", "insert_size_metrics", "flagstats"]),
        opt = expand("{DIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{qc_type}", DIR = CFG["dirs"]["outputs"], qc_type = md.switch_on_wildcard("seq_type", {"_default" : "", "mrna" : "collect_rnaseq_metrics", "genome" : "collect_wgs_metrics", "capture": ["hs_metrics", "interval_hs_metrics"]}))
    output:
        touch(CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{sample_id}.dispatched")


rule _picard_qc_all:
    input:
        vcfs = expand(expand("{dir}{seq_type}--{genome_build}/{sample_id}.dispatched", zip,
                    dir = CFG["dirs"]["outputs"],
                    seq_type = CFG["runs"]["tumour_seq_type"],
                    genome_build = CFG["runs"]["tumour_genome_build"],
                    sample_id = CFG["samples"]["sample_id"]),
                var_type = ['indels', 'snvs'])

##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
