#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["varscan"]`
CFG = op.setup_module(
    name = "varscan",
    version = "1.0",
    subdirectories = ["inputs", "mpileup", "varscan", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _varscan_input_bam,
    _varscan_output_vcf,
    _varscan_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _varscan_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bam + ".bai", output.bam + ".bai")


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _varscan_bam2mpu:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.mpileup"
    log:
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam2mpu.stderr.log"
    params:
        opts = CFG["options"]["mpileup"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["mpileup"]
    resources:
        mem_mb = CFG["mem_mb"]["mpileup"]
    shell:
        op.as_one_line("""
        samtools mpileup {params.opts}
        -f {input.fasta} {input.bam}
        > {output.mpu}
        2> {log.stderr}
        """)


rule _varscan_somatic:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.mpileup",
        normal_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.mpileup"
    output:
        snp = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snp.vcf",
        indel = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/indel.vcf"
    wildcard_constraints:
        pair_status = "matched|unmatched"
    log:
        stdout = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_somatic.stdout.log",
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_somatic.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["somatic"])
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["somatic"]
    resources:
        mem_mb = CFG["mem_mb"]["somatic"]
    shell:
        op.as_one_line("""
        varscan somatic 
        {input.normal_mpu} {input.tumour_mpu} 
        --output-snp {output.snp} --output-indel {output.indel}
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


rule _varscan_unpaired:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.mpileup",
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
    wildcard_constraints:
        pair_status = "no_normal"
    log:
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_{vcf_name}.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["unpaired"]),
        cns = op.switch_on_wildcard("seq_type", {"cns": CFG["options"]["unpaired"]["cns"], "indel": "", "snp": ""})
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["unpaired"]
    resources:
        mem_mb = CFG["mem_mb"]["unpaired"]
    shell:
        op.as_one_line("""
        varscan mpileup2{wildcards.vcf_name} 
        {input.tumour_mpu} 
        {params.opts}
        {params.cns}
        > {output.vcf} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _varscan_output_vcf:
    input:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


def _varscan_get_output(wildcards):
    if wildcards.pair_status == "no_normal":
        return expand(rules._varscan_output_vcf.output.vcf, vcf_name = ["indel", "snp", "cns"], **wildcards)
    else: 
        return expand(rules._varscan_output_vcf.output.vcf, vcf_name = ["indel", "snp"], **wildcards)

rule _varscan_dispatch:
    input:
        _varscan_get_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _varscan_all:
    input:
        expand(
            [
                rules._varscan_dispatch.output.dispatched
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
