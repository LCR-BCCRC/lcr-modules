#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Nicole Thomas
# Module Author:    Bruno Grande
# Contributors:     N/A


##### SETUP #####


# Import standard modules
import os

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["star"]`
CFG = op.setup_module(
    name = "star",
    version = "1.0",
    subdirectories = ["inputs", "star", "sort_bam", "mark_dups", "outputs"],
)

# Include `utils` module
# include: "../../utils/1.0/utils.smk"

# Define rules to be run locally when using a compute cluster
localrules:
    _star_input_fastq,
    _star_symlink_star_bam,
    _star_symlink_sorted_bam,
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
        op.relative_symlink(input.fastq_1, output.fastq_1)
        op.relative_symlink(input.fastq_2, output.fastq_2)


# Align reads using STAR (including soft-clipped chimeric reads)
rule _star_run:
    input:
        fastq_1 = rules._star_input_fastq.output.fastq_1,
        fastq_2 = rules._star_input_fastq.output.fastq_2,
        index = reference_files("genomes/{{genome_build}}/star_index/star-2.7.3a/gencode-{}/overhang-{}".format(
            CFG["reference_params"]["gencode_release"], CFG["reference_params"]["star_overhang"]
        )),
        gtf = reference_files("genomes/{{genome_build}}/annotations/gencode_annotation-{}.gtf".format(
            CFG["reference_params"]["gencode_release"]
        ))
    output:
        bam = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/Aligned.out.bam"
    log:
        stdout = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stdout.log",
        stderr = CFG["logs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/star.stderr.log"
    params:
        opts = CFG["options"]["star"],
        prefix = CFG["dirs"]["star"] + "{seq_type}--{genome_build}/{sample_id}/",
        star_overhang = CFG["reference_params"]["star_overhang"]
    conda:
        CFG["conda_envs"]["star"]
    threads:
        CFG["threads"]["star"]
    resources:
        mem_mb = CFG["mem_mb"]["star"]
    shell:
        op.as_one_line("""
        STAR {params.opts} --readFilesIn {input.fastq_1} {input.fastq_2} --genomeDir {input.index}
        --outFileNamePrefix {params.prefix} --runThreadN {threads} --sjdbGTFfile {input.gtf}
        --sjdbOverhang {params.star_overhang} > {log.stdout} 2> {log.stderr}
            &&
        rmdir {params.prefix}/_STARtmp
        """)


# Create symlink in subdirectory where BAM files will be sorted by the `utils` module
rule _star_symlink_star_bam:
    input:
        bam = rules._star_run.output.bam
    output:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


# Create symlink in subdirectory where duplicates will be marked by the `utils` module
# This rule will trigger the `utils` rule for sorting the BAM file
# By this point, the sorted BAM file exists, so this rule deletes the original BAM file
rule _star_symlink_sorted_bam:
    input:
        bam = CFG["dirs"]["sort_bam"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam",
        star_bam = rules._star_run.output.bam
    output:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        os.remove(input.star_bam)
        shell("touch {input.star_bam}.deleted")


# Symlinks the final output files into the module results directory (under '99-outputs/')
# This rule will trigger the `utils` rules for marking duplicates and indexing the BAM file
# By this point, the mdups BAM file exists, so this rule deletes the sorted BAM file
rule _star_output_bam:
    input:
        bam = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.mdups.bam",
        bai = CFG["dirs"]["mark_dups"] + "{seq_type}--{genome_build}/{sample_id}.sort.mdups.bam.bai",
        sorted_bam = rules._star_symlink_sorted_bam.input.bam
    output:
        bam = CFG["dirs"]["outputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam, in_module = True)
        op.relative_symlink(input.bai, output.bam + ".bai", in_module = True)
        os.remove(input.sorted_bam)
        shell("touch {input.sorted_bam}.deleted")


# Generates the target sentinels for each run, which generate the symlinks
rule _star_all:
    input:
        expand(
            rules._star_output_bam.output.bam,
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
