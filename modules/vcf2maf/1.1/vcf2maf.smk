#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


##### SETUP #####

import sys, os
from os.path import join

import oncopipe as op

# Setup module and store module-specific configuration in `CONFIG`
CFG = op.setup_module(
    name = "vcf2maf",
    version = "1.1",
    subdirectories = ["inputs","decompressed","vcf2maf","crossmap","outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _vcf2maf_input_vcf,
    _vcf2maf_decompress_vcf,
    _vcf2maf_output_maf,
    _vcf2maf_crossmap,
    _vcf2maf_all

VERSION_MAP = {
    "grch37": "GRCh37",
    "hg38": "GRCh38",
    "hs37d5": "GRCh37"
}

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _vcf2maf_input_vcf:
    input:
        vcf_gz = CFG["inputs"]["sample_vcf_gz"]
    output:
        vcf_gz = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz"
    run:
        op.relative_symlink(input.vcf_gz, output.vcf_gz)

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
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf", 
        vep = temp(CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vep.vcf")
    log:
        stdout = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stderr.log",
    params:
        opts = CFG["options"]["vcf2maf"],
        build = lambda w: VERSION_MAP[w.genome_build]
    conda:
        CFG["conda_envs"]["vcf2maf"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        mem_mb = CFG["mem_mb"]["vcf2maf"]
    shell:
        op.as_one_line("""
        if [[ -e {output.maf} ]]; then rm -f {output.maf}; fi;
        if [[ -e {output.vep} ]]; then rm -f {output.vep}; fi;
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


def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

rule _vcf2maf_crossmap:
    input:
        maf = rules._vcf2maf_run.output.maf,
        convert_coord = CFG["inputs"]["convert_coord"],
        chains = get_chain
    output:
        maf =  CFG["dirs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.converted"
    log:
        stdout = CFG["logs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stdout.log",
        stderr = CFG["logs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stderr.log"
    conda:
        CFG["conda_envs"]["crossmap"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        mem_mb = CFG["mem_mb"]["vcf2maf"]
    params:
        out_name = CFG["dirs"]["crossmap"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.converted_",
        chain = lambda w: "hg38ToHg19" if "38" in str({w.genome_build}) else "hg19ToHg38",
        file = ".maf"
    shell:
        op.as_one_line("""
        {input.convert_coord}
        {input.maf}
        {input.chains}
        {params.out_name}{params.chain}{params.file}
        crossmap
        > {log.stdout} 2> {log.stderr}
        && touch {output.maf}
        """)


rule _vcf2maf_output_maf:
    input:
        maf = str(rules._vcf2maf_run.output.maf),
        maf_converted = str(rules._vcf2maf_crossmap.output.maf)
    output:
        maf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_{base_name}.maf"
    params:
        chain = lambda w: "hg38ToHg19" if "38" in str({w.genome_build}) else "hg19ToHg38"
    run:
        op.relative_symlink(input.maf, output.maf)
        op.relative_symlink((input.maf_converted+str("_")+str(params.chain)+str(".maf")), (output.maf[:-4]+str(".converted_")+str(params.chain)+str(".maf")))

# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(str(rules._vcf2maf_output_maf.output.maf), zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
            base_name = [CFG["vcf_base_name"]] * len(CFG["runs"]["tumour_sample_id"]))

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)