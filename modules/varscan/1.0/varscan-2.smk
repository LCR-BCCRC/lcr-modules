#!/usr/bin/env snakemake
# -----------------------------------------------------------------------------
# Author:        Helena Winata
# email/github:  hwinata@bccrc.ca / whelena
# -----------------------------------------------------------------------------
# Input:         {sample_id}.bam  
#
# Output:        {sample_id}/snvs.vcf      
#                {sample_id}/indels.vcf        
#
# Purpose: Processed bam files are converted into mpileup files using samtools.
#          Varscan then call variants (snvs and indels) from mpileup files and store them as VCF files
# Modes: 
#   unpaired: run varscan on unpaired tumour data
#   paired: run varscan on paired tumour-normal or tumour-unmatched normal data
#  -----------------------------------------------------------------------------

import os
import gzip
import modutils as md

##### SETUP #####

CFG = md.setup_module(
    config = config, 
    name = "varscan", 
    version = "1.0",
    subdirs = ["inputs", "mpileup", "vcf", "maf", "outputs"],
    req_references = ["genome_fasta", "vep_dir"],
    scratch_subdirs = ["mpileup"]
)

localrules: 
    _varscan_input, 
    _varscan_output, 
    _varscan_all

include: "vcf2maf.smk"

##### RULES #####

rule _varscan_input:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)
        md.symlink(input.bam + ".bai", output.bam + ".bai")


rule _varscan_bam2mpu:
    input:
        bam = rules._varscan_input.output.bam
    output:
        mpu = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.mpileup")
    log:
        CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam2mpu.stderr.log"
    params:
        opts = CFG["options"]["bam2mpu"],
        fasta  = md.get_reference(CFG, "genome_fasta")
    conda:
        CFG["conda_envs"].get("samtools") or "envs/samtools-0.1.18.yaml"
    threads:
        CFG["threads"].get("bam2mpu") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("bam2mpu") or 8000
    shell:
        md.as_one_line("""
        samtools mpileup {params.opts}
        -f {params.fasta} {input.bam}
        > {output.mpu} 2> {log}
        """)
        

rule _varscan_mpu2vcf_somatic:
    input:
        normalMPU = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.mpileup",
        tumourMPU = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.mpileup"
    output:
        snp = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/snp.vcf",
        indel = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/indel.vcf"
    wildcard_constraints:
        pair_status = "matched|unmatched"
    log:
        stdout = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_somatic.stdout.log",
        stderr = CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_somatic.stderr.log"
    params:
        opts = md.switch_on_wildcard("seq_type", CFG["options"]["somatic"])
    conda:
        CFG["conda_envs"].get("varscan") or "envs/varscan-2.4.4.yaml"
    threads:
        CFG["threads"].get("somatic") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("somatic") or 5000
    shell:
        md.as_one_line("""
        varscan somatic 
        {input.normalMPU} {input.tumourMPU} 
        --output-snp {output.snp} --output-indel {output.indel}
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


rule _varscan_mpu2vcf_unpaired:
    input:
        tumourMPU = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.mpileup"
    output:
        vcf = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
    wildcard_constraints:
        normal_id = "None",
        pair_status = "no_normal"
    log:
        CFG["logs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_{vcf_name}.stderr.log"
    params:
        opts = md.switch_on_wildcard("seq_type", CFG["options"]["unpaired"]),
        cns = md.switch_on_wildcard("vcf_name", {"cns": CFG["options"]["unpaired"]["cns"], "indel": "", "snp": ""})
    conda:
        CFG["conda_envs"].get("varscan") or "envs/varscan-2.4.4.yaml"
    threads:
        CFG["threads"].get("unpaired") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("unpaired") or 5000
    shell:
        md.as_one_line("""
        varscan mpileup2{wildcards.vcf_name} {input.tumourMPU} 
        {params.opts}
        {params.cns}
        > {output.vcf} 2> {log}
        """)


rule _varscan_output:
    input: 
        vcf = CFG["dirs"]["vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf",
        maf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.maf"
    output:
        vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.vcf",
        maf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{vcf_name}/{tumour_id}--{normal_id}--{pair_status}.{vcf_name}.maf"
    run:
        md.symlink(input.vcf, output.vcf) 
        md.symlink(input.maf, output.maf) 

single_vcfs = CFG["outputs"]["unpaired_vcfs"]
def _get_varscan_files(wildcards):
    if wildcards.pair_status == "no_normal":
        vcf_names = single_vcfs
    else:
        vcf_names = ["snp", "indel"]
    vcf_targets = expand(rules._varscan_output.output.vcf, vcf_name = vcf_names, **wildcards)
    return vcf_targets


rule _varscan_all_dispatch:
    input:
        _get_varscan_files
    output:
        touch(CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


rule _varscan_all:
    input:
        expand(rules._varscan_all_dispatch.output, zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"])

'''
        expand(expand("{{dir}}{seq_type}--{genome_build}/{{vcf_name}}/{tumour_id}--{normal_id}--{pair_status}.{{vcf_name}}.vcf", zip,
                    seq_type=CFG["runs"]["tumour_seq_type"],
                    genome_build=CFG["runs"]["tumour_genome_build"],
                    tumour_id=CFG["runs"]["tumour_sample_id"],
                    normal_id=CFG["runs"]["normal_sample_id"],
                    pair_status=CFG["runs"]["pair_status"]),
                vcf_name=["indel", "snp", "cns"],
                dir = CFG["dirs"]["outputs"])
'''

##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
