#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mutect2"]`
CFG = op.setup_module(
    name = "mutect2",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "mutect2", "filter", "decompress", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _mutect2_input_bam,
    _mutect2_output_vcf,
    _mutect2_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutect2_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


# Launces Mutect2 in paired mode
rule _mutect2_run:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.vcf.gz"
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_run.stderr.log"
    params:
        opts = CFG["options"]["mutect2_run"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    resources:
        mem_mb = CFG["mem_mb"]["mutect2_run"]
    shell:
        op.as_one_line("""
        gatk Mutect2 {params.opts} -I {input.tumour_bam} -I {input.normal_bam}
        -R {input.fasta} -normal {wildcards.normal_id} -O {output.vcf} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Marks variants filtered or PASS annotations
rule _mutect2_filter:
    input:
        vcf = str(rules._mutect2_run.output.vcf),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.vcf.gz"
    log:
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_filter.stderr.log"
    params:
        opts = CFG["options"]["mutect2_filter"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_filter"]
    resources:
        mem_mb = CFG["mem_mb"]["mutect2_filter"]
    shell:
        op.as_one_line("""
        gatk FilterMutectCalls {params.opts} -V {input.vcf} -R {input.fasta}
        -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)

# Decompress Mutect2 for downstream analyses
rule _mutect2_decompress:
    input:
        vcf = str(rules._mutect2_filter.output.vcf)
    output:
        vcf = CFG["dirs"]["decompress"] + "{seq_type}--{genome_build}/{tumour_id}--{normal}--{pair_status}/output.filt.vcf"
    log:
        stderr = CFG["logs"]["decompress"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/mutect2_decompress.stderr.log"
    shell:
        "bgzip -d -c {input.vcf} > {output.vcf} 2> {log.stderr}"

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutect2_output_vcf:
    input:
        vcf = str(rules._mutect2_decompress.output.vcf)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutect2_all:
    input:
        expand(
            [
                str(rules._mutect2_output_vcf.output.vcf),
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
