#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Jasper
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import os

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["shimmer"]`
CFG = op.setup_module(
    name = "shimmer",
    version = "1.0",
    subdirectories = ["inputs", "shimmer", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _shimmer_input_bam,
    _shimmer_output_vcf,
    _shimmer_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _shimmer_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


rule _shimmer_on:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_diffs.vcf",
        vs = CFG["dirs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_diffs.vs",
        indels = CFG["dirs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_indels.vs"
    log:
        stdout = CFG["logs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/shimmer.stdout.log",
        stderr = CFG["logs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/shimmer.stderr.log"
    params:
        shimmer = CFG["software"]["shimmer"],
        printCompCounts = CFG["software"]["printCompCounts"],
        dirOut = CFG["dirs"]["shimmer"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/",
        fdr_cutoff = CFG["options"]["fdrcutoff"],
        minBaseQual = CFG["options"]["minBaseQualityScore"],
        tumour_bam = os.path.realpath(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam"),
        normal_bam = os.path.realpath(CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam")
    conda:
        CFG["conda_envs"]["shimmer"]
    threads:
        CFG["threads"]["shimmer"]
    resources:
        mem_mb = CFG["mem_mb"]["shimmer"]
    shell:
        """
        absPath=$(pwd) ;
        export PATH=$absPath/{params.printCompCounts}:$PATH ;
        mkdir -p {params.dirOut} && {params.shimmer} {params.normal_bam} {params.tumour_bam} --ref {input.fasta} --max_q {params.fdr_cutoff} --minqual {params.minBaseQual} --outdir {params.dirOut}  > {log.stdout} 2> {log.stderr}
        """
    


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _shimmer_output_vcf:
    input:
        vcf = str(rules._shimmer_on.output.vcf),
        vs = str(rules._shimmer_on.output.vs),
        indels = str(rules._shimmer_on.output.indels)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.somatic_diffs.vcf",
        vs = CFG["dirs"]["outputs"] + "vs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.somatic_diffs.vs",
        indels = CFG["dirs"]["outputs"] + "indels/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.somatic_indels.vs"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.vs, output.vs)
        op.relative_symlink(input.indels, output.indels)


# Generates the target sentinels for each run, which generate the symlinks
rule _shimmer_all:
    input:
        expand(
            [
                str(rules._shimmer_output_vcf.output.vcf),
                str(rules._shimmer_output_vcf.output.vs),
                str(rules._shimmer_output_vcf.output.indels)
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
