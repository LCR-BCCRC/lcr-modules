#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Cancer IT
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cluster_sv"]`
CFG = op.setup_module(
    name = "cluster_sv",
    version = "1.0",
    subdirectories = ["inputs", "annotate", "combine", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_input_vcf,
    _cluster_sv_decompress_vcf,
    _cluster_sv_filter_vcf,
    _cluster_sv_annotate,
     _cluster_sv_get_map,
    _cluster_sv_combine_bedpe,
    _cluster_sv_all

# Set cohort name to "ALL" if not defined in config file
if not CFG["options"]["cohort"]:
    CFG["options"]["cohort"] = "ALL"

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cluster_sv_input_vcf:
    input:
        vcf_gz = CFG["inputs"]["sample_vcf"]
    output:
        vcf_gz = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf.gz"
    run:
        op.relative_symlink(input.vcf_gz, output.vcf_gz)


rule _cluster_sv_decompress_vcf:
    input:
        vcf_gz = str(rules._cluster_sv_input_vcf.output.vcf_gz)
    output:
        vcf = temp(CFG["dirs"]["annotate"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf")
    shell:
        "gzip -dc {input.vcf_gz} > {output.vcf}" #this should work on both gzip and bcftools compressed files 


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _cluster_sv_filter_vcf:
    input:
        vcf = str(rules._cluster_sv_decompress_vcf.output.vcf)
    output:
        vcf = temp(CFG["dirs"]["annotate"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.vcf"),
    params:
        filter_vcf= CFG["options"]["filter_vcf"]
    conda: 
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        sed -e 's/chr//g' {input.vcf} |
        bedtools intersect -v -header
        -a stdin
        -b {params.filter_vcf}
        > {output.vcf}
        """)

def _cluster_sv_get_vcf(wildcards):
    CFG = config["lcr-modules"]["cluster_sv"]
    RUNS = config["lcr-modules"]["cluster_sv"]["runs"]

    if CFG["options"]["filter_vcf"] is None:
        print('No filter is applied to VCF file')
        vcf_file = expand(
            [
                str(rules._cluster_sv_input_vcf.output.vcf),
            ],
            zip,
            **wildcards)
    else:
        vcf_file = expand(
            [
                str(rules._cluster_sv_filter_vcf.output.vcf),
            ],
            zip,
            **wildcards)
    
    return { 'vcf': vcf_file }


rule _cluster_sv_annotate:
    input:
        unpack(_cluster_sv_get_vcf)
    output:
        vcf = CFG["dirs"]["annotate"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.vcf",
        bedpe = CFG["dirs"]["annotate"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.bedpe"
    log: 
        CFG["logs"]["annotate"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_R_annotateSV.log"
    params:
        r_script = CFG["options"]["vcf_annotation_script"]
    conda: 
        CFG["conda_envs"]["gridss"]
    shell:
        op.as_one_line("""
        Rscript {params.r_script} 
        -vcf {input.vcf}
        -ref {wildcards.genome_build}
        -out_vcf {output.vcf} 
        -out_bedpe {output.bedpe} 
        &> {log}
        """)

'''
rule _cluster_sv_vcf2bedpe:
    input:
        vcf = str(rules._cluster_sv_annotate.output.vcf)
    output:
        bedpe = CFG["dirs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.bedpe",
    log: 
        CFG["logs"]["bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.log"
    conda: 
        CFG["conda_envs"]["svtools"]
    shell:
        op.as_one_line("""
        svtools vcftobedpe -i {input.vcf} -o {output.bedpe}
        """)



rule _cluster_sv_reformat_bedpe:
    input:
        bedpe = str(rules._cluster_sv_annotate.output.bedpe)
    output:
        bedpe = temp(CFG["dirs"]["reformat_bedpe"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"),
        mapping = temp(CFG["dirs"]["reformat_bedpe"] + "mapping/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_map.tsv")
    shell:
        op.as_one_line("""
        sed '/^#/d' {input.bedpe} |
        sed -e 's/chr//g' > "{output.bedpe}"
        &&
        awk -F '\t' 'BEGIN {{OFS="\t"}};
        {{ print "{wildcards.tumour_id}", "{wildcards.normal_id}", $7 }}'
        > {output.mapping}
        """)
'''

rule _cluster_sv_get_map:
    input:
        bedpe = str(rules._cluster_sv_annotate.output.bedpe)
    output:
        #bedpe = temp(CFG["dirs"]["reformat_bedpe"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"),
        mapping = temp(CFG["dirs"]["annotate"] + "mapping/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_map.tsv")
    shell:
        op.as_one_line("""
        awk -F '\t' 'BEGIN {{OFS="\t"}};
        {{ print "{wildcards.tumour_id}", "{wildcards.normal_id}", $7 }}'
        {input.bedpe} > {output.mapping}
        """)

def _cluster_sv_get_bedpe(wildcards):
    CFG = config["lcr-modules"]["cluster_sv"]
    SAMPLES = config["lcr-modules"]["cluster_sv"]["samples"]
    RUNS = config["lcr-modules"]["cluster_sv"]["runs"]

    # if wildcards.cohort != 'ALL':
    if wildcards.cohort in SAMPLES.cohort:
        TUMOURS = SAMPLES.loc[SAMPLES.cohort == wildcards.cohort].sample_id
        RUNS = RUNS.loc[RUNS.tumour_sample_id.isin(TUMOURS)]

    bedpe_files = expand(
        [
            str(rules._cluster_sv_annotate.output.bedpe),
        ],
        zip,
        seq_type=RUNS["tumour_seq_type"],
        genome_build=RUNS["tumour_genome_build"],
        tumour_id=RUNS["tumour_sample_id"],
        normal_id=RUNS["normal_sample_id"],
        pair_status=RUNS["pair_status"])

    map_files = expand(
        [
            str(rules._cluster_sv_get_map.output.mapping),
        ],
        zip, 
        seq_type=RUNS["tumour_seq_type"],
        genome_build=RUNS["tumour_genome_build"],
        tumour_id=RUNS["tumour_sample_id"],
        normal_id=RUNS["normal_sample_id"],
        pair_status=RUNS["pair_status"])

    return { 'bedpe': bedpe_files, 'mapping': map_files }


rule _cluster_sv_combine_bedpe:
    input:
        unpack(_cluster_sv_get_bedpe)
    output:
        bedpe = CFG["dirs"]["combine"] + "combined/{seq_type}--{genome_build}/{cohort}_combined.bedpe",
        mapping = CFG["dirs"]["combine"] + "combined/{seq_type}--{genome_build}/{cohort}_map.tsv"
    shell:
        op.as_one_line("""
        cat {input.bedpe} > {output.bedpe} &&
        cat {input.mapping} > {output.mapping} &&
        sed -i '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tALT_A\tALT_B' {output.bedpe} &&
        sed -i '1i #tumour_id\tnormal_id\tSV_id' {output.mapping}
        """)


rule _cluster_sv_all:
    input:
        expand(
            [
                str(rules._cluster_sv_combine_bedpe.output.bedpe),
                str(rules._cluster_sv_combine_bedpe.output.mapping)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            cohort=CFG["options"]["cohort"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
