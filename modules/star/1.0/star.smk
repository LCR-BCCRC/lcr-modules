#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Snakefile Author:    Nicole Thomas
# Module Author:                Bruno Grande
# Additional Contributors:      N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["star"]`
CFG = op.setup_module(
    name = "star",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "star", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _star_input_fastq,
    _star_step_2,
    _star_output_bam,
    _star_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _star_input_fastq:
    input:
        fastq = CFG["inputs"]["sample_fastq"]
    output:
        fastq = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.fastq"
    run:
        op.symlink(input.fastq, output.fastq)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _star_step_1:
    input:
        fastq = rules._star_input_fastq.output.fastq,
        fasta = op.get_reference(CFG, "genome_fasta")
    output:
        bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/output.bam"
    log:
        stdout = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["star"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --input {input.fastq} --ref-fasta {params.fasta}
        --output {output.bam} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _star_step_2:
    input:
        bam = rules._star_step_1.output.bam
    output:
        bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/output.filt.bam"
    log:
        stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.bam} > {output.bam} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _star_output_bam:
    input:
        bam = rules._star_step_2.output.bam
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.output.filt.bam"
    run:
        op.symlink(input.bam, output.bam)


# Generates the target sentinels for each run, which generate the symlinks
rule _star_all:
    input:
        expand(
            [
                rules._star_output_bam.output.bam,
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
