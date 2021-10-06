#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["ega_download"]`
CFG = op.setup_module(
    name = "ega_download",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "ega_download", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _ega_download_input_csv,
    _ega_download_step_2,
    _ega_download_output_bam,
    _ega_download_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _ega_download_input_csv:
    input:
        csv = CFG["inputs"]["sample_csv"]
    output:
        csv = CFG["dirs"]["inputs"] + "csv/{seq_type}--{genome_build}/{sample_id}.csv"
    run:
        op.relative_symlink(input.csv, output.csv)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _ega_download_step_1:
    input:
        csv = str(rules._ega_download_input_csv.output.csv),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = CFG["dirs"]["ega_download"] + "{seq_type}--{genome_build}/{sample_id}/output.bam"
    log:
        stdout = CFG["logs"]["ega_download"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["ega_download"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
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
        <TODO> {params.opts} --input {input.csv} --ref-fasta {input.fasta}
        --output {output.bam} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _ega_download_step_2:
    input:
        bam = str(rules._ega_download_step_1.output.bam)
    output:
        bam = CFG["dirs"]["ega_download"] + "{seq_type}--{genome_build}/{sample_id}/output.filt.bam"
    log:
        stderr = CFG["logs"]["ega_download"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.bam} > {output.bam} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _ega_download_output_bam:
    input:
        bam = str(rules._ega_download_step_2.output.bam)
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.output.filt.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


# Generates the target sentinels for each run, which generate the symlinks
rule _ega_download_all:
    input:
        expand(
            [
                str(rules._ega_download_output_bam.output.bam),
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
