#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "ginkgo", 
    version = "1.0",
    # TODO: If applicable, add other output subdirectories
    subdirs = ["inputs", "ginkgo", "outputs"],
    # TODO: Replace "genome_fasta" with actual reference requirements
    req_references = ["genome_fasta"] 
)

# Define rules to be run locally when using a compute cluster.
# TODO: Replace with actual rule names
localrules: 
    _ginkgo_input_bed,
    _ginkgo_output,
    _ginkgo_all


##### RULES #####
# in single cell modules, sample_id and lib_id refers to patient_id and sample_id respectively on the datasheet
# for ginkgo purposes, each library = a cell

rule _ginkgo_input_bed:
    input:
        bed = CFG["inputs"]["sample_bed"]
    output:
        bed = CFG["dirs"]["inputs"] + "{genome_build}/{sample_id}/{lib_id}.bed.gz"
    run:
        md.symlink(input.bed, output.bed)


rule _ginkgo_link_bed_to_bins:
    input:
        bed = rules._ginkgo_input_bed
    output:
        bed = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/{lib_id}.bed.gz"
    run:
        md.symlink(input.bed, output.bed)


def get_bed_files(wildcards):
    sub_df = md.filter_samples(CFG["samples", patient_id = wildcards.sample_id])
    libs = list(sub_df["sample_id"])
    bed = expand("{dir}{{genome_build}}_bin{{bin}}/{sample_id}/{lib_id}.bed.gz", dir = CFG["dirs"]["ginkgo"], sample_id = wildcards.sample_id, lib_id = libs)
    return bed


rule _ginkgo_run:
    input:
        bed = get_bed_files
    output:
        config["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample}.done"
    log: 
        stdout = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stderr.log"
    params:
        ginkgo = CFG["software"],
        bedDir = CFG["dirs"]["ginkgo"] + "_{bin}/{sample}",
        genome = CFG["options"]["genome"],
        binning = CFG["options"]["binMeth"],
        clustDist = CFG["options"]["distMeth"],
        clustLink = CFG["options"]["clustMeth"],
        opts = CFG["options"]["opts"]
    threads: 2
    resources:
        mem_mb = 8000
    shell:
        as_one_line("""
        echo 
        {params.ginkgo} 
        --input {params.bedDir} 
        --genome {params.genome} 
        --binning {params.binning}
        {params.clustDist} {params.clustLink}
        {params.opts} 
        &> {log}
        && touch {output}
        """)



    output:
        vcf = CFG["dirs"]["ginkgo"] + "{seq_type}--{genome_build}/{sample_id}/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["ginkgo"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["ginkgo"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
        fasta = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("ginkgo") or "envs/ginkgo.yaml"
    threads:
        CFG["threads"]["step_1"]
    resources: 
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        md.as_one_line("""
        <TODO> {params.opts} --bam {input.bam} --ref-fasta {params.fasta} 
        --output {output.vcf} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _ginkgo_step_2:
    input:
        vcf = rules._ginkgo_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["ginkgo"] + "{seq_type}--{genome_build}/{sample_id}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["ginkgo"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/').
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _ginkgo_output_vcf:
    input:
        vcf = rules._ginkgo_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.variants.filt.vcf"
    run:
        md.symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_ginkgo_output_*` rule
rule _ginkgo_all:
    input:
        expand(rules._ginkgo_output.output.vcf, 
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["samples"]["seq_type"],
               genome_build=CFG["samples"]["genome_build"],
               sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
