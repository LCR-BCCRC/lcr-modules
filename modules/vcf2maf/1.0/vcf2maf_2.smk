#!/usr/bin/env snakemake

##### MODULES #####

import os
import modutils as md


##### SETUP #####

VCF = md.setup_module(
    config = config, 
    name = "vcf2maf", 
    version = "1.0",
    subdirs = ["vcf2maf"],
    req_references = ["genome_fasta"]
)

localrules: vcf2maf_all


##### RULES #####
"""
rule vcf2maf_input:
    input:
        vcf = VCF["inputs"].get("vcf")
    output:
        vcf = join(VCF["dirs"]["inputs"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.vcf")
    run:
        symlink(input.vcf, output.vcf)
"""

rule vcf2maf_run:
    input:
        vcf = "{vcf_path}/{vcf_name}.vcf"
    output:
        maf = join(VCF["dirs"]["vcf2maf"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.maf")
    log:
        join(VCF["dirs"]["vcf2maf"], "{seq_type}", 
             "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.log.txt")
    params:
        opts = VCF["options"]["vcf2maf"],
        fasta  = config["reference"]["genome_fasta"],
        vep = VCF["inputs"]["vep"]
    conda:
        VCF["conda_envs"]["vcf2maf"] or "envs/vcf2maf.yaml"
    threads:
        VCF["threads"].get("vcf2maf") or 1
    resources: 
        mem_mb = VCF["memory"].get("vcf2maf") or 4000
    shell:
        """
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --ref-fasta {params.fasta} --vep-data {params.vep} --vep-path $vepPATH {params.opts} 2> {log}
        """
        

rule vcf2maf_output:
    input:
        maf = rules.vcf2maf_run.output.maf
    output:
        maf = join(VCF["dirs"]["outputs"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.maf")
    run:
        symlink(input.maf, output.maf)

rule _vcf2_all_dispatch:
    input:

    output:
        VCF["dirs"]["outputs"] + "/{seq_type}/{tumour_id}--{normal_id}--{is_matched}/{caller}.dispatched"

rule _vcf2maf_all:
    input:
        expand(expand(rules._vcf2_all_dispatch, zip,
                    seq_type = VCF["runs"]["tumour_seq_type"],
                    tumour_id = VCF["runs"]["tumour_sample_id"],
                    normal_id = VCF["runs"]["normal_sample_id"],
                    is_matched = VCF["runs"]["is_matched"]),
                caller = VCFcaller)

rule vcf2maf_all:
    input:
        vcfs = expand(expand("{dir}/{seq_type}/{tumour_id}--{normal_id}--{is_matched}/{{caller}}/{{var_type}}.maf", zip,
                    dir = VCF["dirs"]["outputs"],
                    seq_type = VCF["runs"]["tumour_seq_type"],
                    tumour_id = VCF["runs"]["tumour_sample_id"],
                    normal_id = VCF["runs"]["normal_sample_id"],
                    is_matched = VCF["runs"]["is_matched"]),
                var_type = ['indels', 'snvs'], caller = VCFcaller)


##### CLEANUP #####

cleanup_module(VCF)

del VCF
