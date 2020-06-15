#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["battenberg"]`
CFG = op.setup_module(
    name = "battenberg",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "battenberg", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _battenberg_input_bam,
    _battenberg_step_2,
    _battenberg_output_seg,
    _battenberg_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _battenberg_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _battenberg_step_1:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files(CFG["reference"]["genome_fasta"])
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.seg"
    log:
        stdout = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
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
        <TODO> {params.opts} --tumour {input.tumour_bam} --normal {input.normal_bam}
        --ref-fasta {params.fasta} --output {output.seg} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _battenberg_step_2:
    input:
        seg = rules._battenberg_step_1.output.seg
    output:
        seg = CFG["dirs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.seg"
    log:
        stderr = CFG["logs"]["battenberg"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.seg} > {output.seg} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _battenberg_output_seg:
    input:
        seg = rules._battenberg_step_2.output.seg
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.seg"
    run:
        op.relative_symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
rule _battenberg_all:
    input:
        expand(
            [
                rules._battenberg_output_seg.output.seg,
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
