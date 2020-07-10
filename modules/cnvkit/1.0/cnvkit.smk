#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

import glob

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cnvkit"]`
CFG = op.setup_module(
    name = "cnvkit",
    version = "1.0",
    subdirectories = ["inputs", "batch", "segment", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _cnvkit_input_bam,
    _cnvkit_output_vcf,
    _cnvkit_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cnvkit_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _cnvkit_batch:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        access = reference_files("genomes/{genome_build}/cnvkit_access/cnvkit-0.9.7/access.bed"),
        refFlat = reference_files("genomes/{genome_build}/annotations/refFlat_gencode-33.txt"),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        cnr = CFG["dirs"]["batch"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cnr",
        normal_ref = CFG["dirs"]["batch"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/normal_reference.cnn"
    log:
        stdout = CFG["logs"]["batch"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/batch.stdout.log",
        stderr = CFG["logs"]["batch"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/batch.stderr.log"
    params:
        opts = op.switch_on_wildcards("seq_type", CFG["options"]["batch"])
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["batch"]
    resources:
        mem_mb = CFG["mem_mb"]["batch"]
    shell:
        op.as_one_line("""
        cnvkit.py batch
        {input.tumour_bam}
        --normal {input.normal_bam}
        --access {params.access}
        --annotate {imput.refFlat}
        --fasta {params.fasta}
        --processes {threads}
        --output-reference {output.normal_ref}
        --output-dir $( dirname {output.cnr})
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


rule _cnvkit_segment:
    input:
        cnr = rules._cnvkit_batch.output.cnr
    output:
        cns = CFG["dirs"]["segment"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cns"
    log:
        stdout = CFG["logs"]["segment"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/segment.stdout.log",
        stderr = CFG["logs"]["segment"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/segment.stderr.log"
    params:
        opts = CFG["options"]["segment"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["segment"]
    resources:
        mem_mb = CFG["mem_mb"]["segment"]
    shell:
        op.as_one_line("""
        cnvkit.py segment
        -p {threads}
        {params.opts}
        {input.cnr}
        -o {output.cns}
        > {log.stdout} 2> {log.stderr}
        """)


rule _cnvkit_call:
    input:
        cns = rules._cnvkit_segment.output.cns
    output:
        cns = CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.call.cns"
    log:
        stdout = CFG["logs"]["call"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/call.stdout.log",
        stderr = CFG["logs"]["call"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/call.stderr.log"
    params:
        opts = CFG["options"]["call"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    threads:
        CFG["threads"]["call"]
    resources:
        mem_mb = CFG["mem_mb"]["call"]
    shell:
        op.as_one_line("""
        cnvkit.py call
        -p {threads}
        {params.opts}
        {input.cns}
        -o {output.cns}
        > {log.stdout} 2> {log.stderr}
        """)

# PLOTTING
rule _cnvkit_scatter:
    input:
        cns = rules._cnvkit_segment.output.cns,
        cnr = rules._cnvkit_batch.output.cnr,
        chroms = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        scatter = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_scatter.pdf"
    log:
        stderr = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/scatter.stderr.log"
    params:
        opts = CFG["options"]["scatter"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    shell:
        op.as_one_line("""
        cnvkit.py scatter
        {input.cnr}
        --segment {input.cns}
        --output {output.scatter}
        {params.opts}
        &> {log.stderr}
        """)


def _get_heatmap_samples(wildcards):
    CFG = config["lcr-modules"]["cnvkit"]
    pattern = CFG["dirs"]["call"] + f"{wildcards.seq_type}--{wildcards.genome_build}/*/{wildcards.seq_type}


CFG["dirs"]["call"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.call.cns"

def _get_heatmap_samples(metrics_dir):
    DIR = metrics_dir
    CFG = config["lcr-modules"]["cnvkit"]
    def _get_heatmap_samples_custom(wildcards):
        # retrieve from CFG["samples"] if specifed, otherwise default to shared smaples
        sample = CFG.get("samples") or config["lcr-modules"]["_shared"]["samples"]
        # filter samples by seq_type and genome_build
        fsample = op.filter_samples(sample, seq_type=wildcards.seq_type, genome_build=wildcards.genome_build)
        return expand("{dir}{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.call.cns", 
        dir = DIR, 
        sample_id = list(fsample["sample_id"]), **wildcards)
    return _get_sample_metrics_custom


rule _cnvkit_heatmap:
    input:
        cns = expand(rules._cnvkit_call.output.cns, zip,
            tumour_id=CFG["runs"]["tumour_id"],
            normal_id=CFG["runs"]["normal_id"],
            pair_status=CFG["runs"]["pair_status"])
    output:
        hm = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/merged_heatmap.pdf"
    log:
        stderr = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/heatmap.stderr.log"
    params:
        opts = CFG["options"]["heatmap"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    shell:
        op.as_one_line("""
        cnvkit.py heatmap
        {input.cnr}
        --segment {input.cns}
        --output {output.heatmap}
        {params.opts}
        &> {log.stderr}
        """)


rule _cnvkit_heatmap:
    input:
        cns = rules._cnvkit_segment.output.cns,
        cnr = rules._cnvkit_batch.output.cnr,
        chroms = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        diagram = CFG["dirs"]["plots"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_heatmap.pdf"
    log:
        stderr = CFG["logs"]["plots"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/heatmap.stderr.log"
    params:
        opts = CFG["options"]["heatmap"]
    conda:
        CFG["conda_envs"]["cnvkit"]
    shell:
        op.as_one_line("""
        cnvkit.py heatmap
        {input.cnr}
        --segment {input.cns}
        --output {output.diagram}
        {params.opts}
        &> {log.stderr}
        """)


rule _cnvkit_step_2:
    input:
        vcf = rules._cnvkit_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.vcf"
    log:
        stderr = CFG["logs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.vcf} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _cnvkit_output_vcf:
    input:
        vcf = rules._cnvkit_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
rule _cnvkit_all:
    input:
        expand(
            [
                rules._cnvkit_output_vcf.output.vcf,
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
