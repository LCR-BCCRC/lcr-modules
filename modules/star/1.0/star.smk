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
    subdirectories = ["inputs", "star", "sort_bam", "mark_dups", "outputs"],
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
rule _star_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R2.fastq.gz",
    run:
        op.symlink(input.fastq_1, output.fastq_1)
        op.symlink(input.fastq_2, output.fastq_2)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _star_run:
    input:
        fastq = rules._star_input_fastq.output.fastq,
        fasta = op.get_reference(CFG, "star_index")
    output:
        bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/output.bam"
    log:
        stdout = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stdout.log",
        stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stderr.log"
    params:
        opts = CFG["options"]["star"]
    conda:
        CFG["conda_envs"]["star"]
    threads:
        CFG["threads"]["star"]
    resources:
        mem_mb = CFG["mem_mb"]["star"]
    shell:
        op.as_one_line("""
        STAR --runMode alignReads --runThreadN {threads} --genomeDir {input.STAR_index} --genomeLoad NoSharedMemory
        --sjdbGTFfile {input.gtf} --sjdbOverhang 74 --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat
        --outFileNamePrefix {params.STAR_prefix} --outSAMtype BAM SortedByCoordinate
        --outSAMattrIHstart 0 --chimOutType WithinBAM SoftClip --chimSegmentMin 20 --twopassMode Basic
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
