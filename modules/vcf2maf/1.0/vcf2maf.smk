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
CONFIG = config["lcr-modules"]["vcf2maf"]
PATH = ".*\/(" + "|".join(CONFIG["paired_modules"]) + ").*\/"
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"


##### RULES #####

rule:
    input:
        vcf = "{out_dir}/{prefix}/{suffix}.vcf",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CONFIG["inputs"]["vep_cache"]
    output:
        maf = "{out_dir}/{prefix}/{suffix}.maf",
    log:
        stdout = "{out_dir}" + LOG + "/{prefix}/{suffix}_vcf2maf.stdout.log",
        stderr = "{out_dir}" + LOG + "/{prefix}/{suffix}_vcf2maf.stderr.log",
    params:
        opts = CONFIG["options"]["vcf2maf"]
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
        --tumor-id {wildcards.tumour_id} --normal_id {wildcards.normal_id}
        --vcf-tumor-id TUMOR --vcf-normal_id NORMAL
        --ref-fasta {params.fasta}
        --vep-data {params.vep} --vep-path $vepPATH {params.opts} 2> {log}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _vcf2maf_output_maf:
    input:
        maf = rules._vcf2maf_run.output.maf
    output:
        maf = CONFIG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{caller}/{vartype}.maf"
    run:
        op.relative_symlink(input, output)


def _vcf2maf_get_output(wildcards):
    CONFIG = config["lcr-modules"]["vcf2maf"]



# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(rules._vcf2maf_output_maf.output.maf,
            zip,  # Run expand() with zip(), not product()
            seq_type=CONFIG["runs"]["tumour_seq_type"],
            genome_build=CONFIG["runs"]["tumour_genome_build"],
            tumour_id=CONFIG["runs"]["tumour_sample_id"],
            normal_id=CONFIG["runs"]["normal_sample_id"],
            pair_status=CONFIG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CONFIG` variable
op.cleanup_module(CONFIG)
