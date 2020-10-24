#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["starfish"]`
CFG = op.setup_module(
    name = "starfish",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "starfish", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _starfish_input_vcf,
    _starfish_step_2,
    _starfish_output_vcf,
    _starfish_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _starfish_input_vcf:
    input:
        vcf = CFG["inputs"]["sample_vcf"]
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _starfish_step_1:
    input:
        vcf = str(rules._starfish_input_vcf.output.vcf),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf"
    log:
        stdout = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --input {input.vcf} --ref-fasta {input.fasta}
        --output {output.vcf} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _starfish_step_2:
    input:
        vcf = str(rules._starfish_step_1.output.vcf)
    output:
        vcf = CFG["dirs"]["starfish"] + "{seq_type}--{genome_build}/{sample_id}/output.filt.vcf"
    log:
        stderr = CFG["logs"]["starfish"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.vcf} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _starfish_output_vcf:
    input:
        vcf = str(rules._starfish_step_2.output.vcf)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.output.filt.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
rule _starfish_all:
    input:
        expand(
            [
                str(rules._starfish_output_vcf.output.vcf),
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
