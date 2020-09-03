#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####

import sys, os
from os.path import join

import oncopipe as op

# Setup module and store module-specific configuration in `CONFIG` and use LOG to specify log directory
CONFIG = config["lcr-modules"]["vcf2maf"]
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"

VERSION_UPPER = {
    "grch37": "GRCh37",
    "GRCh37": "GRCh37",
    "grch38": "GRCh38",
    "GRCh38": "GRCh38",
}

localrules: _vcf2maf_decompress_vcf

##### RULES #####
rule _vcf2maf_decompress_vcf:
    input: 
        vcf = "{out_dir}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz"
    output: 
        vcf = temp("{out_dir}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf")
    wildcard_constraints:
        out_dir = ".*(maf).*"
    shell:
        "bgzip -dc {input.vcf} > {output.vcf}"


rule _vcf2maf_run:
    input:
        vcf = "{out_dir}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CONFIG["inputs"]["vep_cache"]
    output:
        maf = "{out_dir}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.maf",
    log:
        stdout = "{out_dir}" + LOG + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}_vcf2maf.stdout.log",
        stderr = "{out_dir}" + LOG + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}_vcf2maf.stderr.log",
    wildcard_constraints:
        out_dir = ".*(maf).*"
    params:
        opts = CONFIG["options"]["vcf2maf"],
        build = lambda w: f"{VERSION_UPPER[w.genome_build]}"
    conda:
        CONFIG["conda_envs"]["vcf2maf"]
    threads:
        CONFIG["threads"]["vcf2maf"]
    resources:
        mem_mb = CONFIG["mem_mb"]["vcf2maf"]
    shell:
        op.as_one_line("""
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*;
        vcf2maf.pl 
        --input-vcf {input.vcf} 
        --output-maf {output.maf} 
        --tumor-id {wildcards.tumour_id} --normal-id {wildcards.normal_id}
        --ref-fasta {input.fasta}
        --ncbi-build {params.build}
        --vep-data {input.vep_cache}
        --vep-path $vepPATH {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)

