#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####

import sys, os
from os.path import join

import oncopipe as op

# Setup module and store module-specific configuration in `CONFIG`
CFG = op.setup_module(
    name = "vcf2maf",
    version = "1.0",
    subdirectories = ["inputs","decompressed","vcf2maf","outputs"]
)


wildcard_constraints:
    prefix = "[0-9]{2}-.*"

VERSION_UPPER = {
    "grch37": "GRCh37",
    "GRCh37": "GRCh37",
    "grch38": "GRCh38",
    "GRCh38": "GRCh38",
}

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _vcf2maf_input_vcf:
    params:
        vcf_gz = CFG["inputs"]["sample_vcf_gz"]
    output:
        vcf_gz = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz",
    run:
        op.relative_symlink(params.vcf_gz, output.vcf_gz)

rule _vcf2maf_decompress_vcf:
    input:
        vcf_gz = str(rules._vcf2maf_input_vcf.output.vcf_gz)
    output:
        vcf = CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf"
    shell:
        "gzip -dc {input.vcf_gz} > {output.vcf}" #this should work on both gzip and bcftools compressed files 

rule _vcf2maf_run:
    input:
        vcf = CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CFG["inputs"]["vep_cache"]
    output:
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf"
    log:
        stdout = CFG["logs"]["vcf2maf"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stderr.log",
    params:
        opts = CFG["options"]["vcf2maf"]
    conda:
        CFG["conda_envs"]["vcf2maf"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        mem_mb = CFG["mem_mb"]["vcf2maf"]
    shell:
        op.as_one_line("""
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*;
        vcf2maf.pl 
        --input-vcf {input.vcf} 
        --output-maf {output.maf} 
        --tumor-id {wildcards.tumour_id} --normal-id {wildcards.normal_id}
        --ref-fasta {input.fasta}
        --ncbi-build {wildcards.genome_build}
        --vep-data $(dirname $(dirname {input.vep_cache}))
        --vep-path $vepPATH {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)

# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(str(rules._vcf2maf_run.output.maf), zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
            base_name = CFG["vcf_base_name"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)