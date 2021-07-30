#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Cancer IT
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####

import sys, os
from os.path import join

# Import package with useful functions for developing analysis modules
import oncopipe as op

CFG = config["lcr-modules"]["cluster_sv"]
LOG = "/logs/" + op._session.launched_fmt

wildcard_constraints:
    prefix = "[0-9]{2}-.*"

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_combine_bed_filter,
    _cluster_sv_bed_filter,
    _cluster_sv_blacklist_filter,
    _cluster_sv_ffpe_filter,
    _cluster_sv_symlink_filter
    _cluster_sv_get_map,
    _cluster_sv_combine_bedpe,
    _cluster_sv_combined_all


##### RULES #####
rule _cluster_sv_combine_bed_filter:
    input:
        bed_list = CFG["options"]["filter"]["bed_list"]
    output:
        bed = CFG["dirs"]["filter"] + "bed_filter/{seq_type}--{genome_build}/combined_filter.bed"
    shell:
        op.as_one_line("""
        cat {input.bed_list} > {output.bed}
        """)

# Filter bedpe file with combined bed file (gnomad and germline from Jabba PON)
rule _cluster_sv_bed_filter:
    input:
        bedpe = str(rules._cluster_sv_annotate.output.bedpe),
        bed = str(rules._cluster_sv_combine_bed_filter.output.bed)
    output:
        bedpe = CFG["dirs"]["filter"] + "bed_filter/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bed_filtered.bedpe"
    conda: 
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        sed -e 's/chr//g' {input.bedpe} |
        pairToBed -a stdin
        -b {input.bed} -bedpe -type neither 
        > {output_bedpe}
        """)


rule _cluster_sv_blacklist_filter:
    input:
        bedpe = str(rules._cluster_sv_bed_filter.output.bedpe)
    output:
        bedpe = CFG["dirs"]["filter"] + "blacklist_filter/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.blacklisted.vcf"
    params:
        blacklist= CFG["options"]["filter"]["blacklist"],
        args = CFG["options"]["blacklist_filter"]
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        pairToPair -a {input.bedpe}
        -b {params.blacklist} 
        -type neither
        {params.args}
        > {output.bedpe}
        """)

rule _cluster_sv_ffpe_filter:
    input:
        bedpe = str(rules._cluster_sv_blacklist_filter.output.bedpe)
    output:
        bedpe = CFG["dirs"]["filter"] + "FFPE/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.ffpe_filtered.bedpe"
    params:
        r_script = CFG["scripts"]["ffpe_filter_script"],
        args = CFG["options"]["ffpe_filter"]
    conda: 
        CFG["conda_envs"]["gridss"]
    shell:
        op.as_one_line("""
        Rscript {params.r_script}
        -bedpe {input.bedpe}
        {params.args}
        -out_bedpe {output.bedpe}
        """)

ffpe_switches = {"_default": str(rules._cluster_sv_blacklist_filter.output.bedpe),
                 "FFPE": str(rules._cluster_sv_ffpe_filter.output.bedpe)}

rule _cluster_sv_symlink_filter:
    input:
        bedpe = op.switch_on_column("ffpe_or_frozen", CFG["samples"], ffpe_switches, match_on = "tumour")
    output:
        bedpe = CFG["dirs"]["filter"] + "final/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)

#######################################################################
# GET COMBINED BEDPE AND MAP
#######################################################################
# Set cohort name to "ALL" if not defined in config file
if not CFG["options"]["cohort"]:
    CFG["options"]["cohort"] = "ALL"

rule _cluster_sv_get_map:
    input:
        bedpe = str(rules._cluster_sv_symlink_filter.output.bedpe)
    output:
        mapping = temp(CFG["dirs"]["combine"] + "mapping/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_map.tsv")
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
            str(rules._cluster_sv_symlink_filter.output.bedpe)
        ],
        zip,
        seq_type=RUNS["tumour_seq_type"],
        genome_build=RUNS["tumour_genome_build"],
        tumour_id=RUNS["tumour_sample_id"],
        normal_id=RUNS["normal_sample_id"],
        pair_status=RUNS["pair_status"])

    map_files = expand(
        [
            str(rules._cluster_sv_get_map.output.mapping)
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
        bedpe = CFG["dirs"]["combine"] + "{seq_type}--{genome_build}/{cohort}_combined.bedpe",
        mapping = CFG["dirs"]["combine"] + "{seq_type}--{genome_build}/{cohort}_map.tsv"
    shell:
        op.as_one_line("""
        cat {input.bedpe} > {output.bedpe} &&
        cat {input.mapping} > {output.mapping} &&
        sed -i '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tALT_A\tALT_B' {output.bedpe} &&
        sed -i '1i #tumour_id\tnormal_id\tSV_id' {output.mapping}
        """)


# get combined bedpe and mapping file
rule _cluster_sv_combined_all:
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

