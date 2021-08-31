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
    subdirectories = ["inputs", "annotate", "filter",
                      "combine", "cluster_sv", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_install,
    _cluster_sv_input_vcf,
    _cluster_sv_decompress_vcf,
    _cluster_sv_annotate,
    # FILTERING RULES
    _cluster_sv_combine_bed_filter,
    _cluster_sv_bed_filter,
    _cluster_sv_blacklist_filter,
    _cluster_sv_ffpe_filter,
    _cluster_sv_symlink_filter,
     # COMBINING RULES
    _cluster_sv_get_map,
    _cluster_sv_combine_bedpe,
    _cluster_sv_combined_all,
    _cluster_sv_output_tsv,
    _cluster_sv_all


##### RULES #####
rule _cluster_sv_install:
    output:
        cluster_sv = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/R/run_cluster_sv.R" # main R script from the repo
    params:
        url = "https://github.com/whelena/ClusterSV/archive/refs/tags/v" + str(CFG["options"]["cluster_sv_version"]) + ".tar.gz",
        folder = CFG["dirs"]["inputs"]
    shell:
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
        """)


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
        vcf = temp(CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf")
    shell:
        "gzip -dc {input.vcf_gz} > {output.vcf}" #this should work on both gzip and bcftools compressed files 


rule _cluster_sv_annotate:
    input:
        vcf = str(rules._cluster_sv_decompress_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["annotate"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.vcf",
        bedpe = CFG["dirs"]["annotate"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.bedpe"
    log: 
        CFG["logs"]["annotate"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_R_annotateSV.log"
    params:
        r_script = CFG["scripts"]["vcf_annotation_script"],
        vaf_thres = CFG["options"]["vaf_threshold"]
    conda: 
        CFG["conda_envs"]["gridss"]
    shell:
        op.as_one_line("""
        Rscript {params.r_script}
        -vcf {input.vcf}
        -ref {wildcards.genome_build}
        -vaf_thres {params.vaf_thres}
        -out_vcf {output.vcf}
        -out_bedpe {output.bedpe}
        &> {log}
        """)


#######################################################################
# FILTER BEDPE
#######################################################################
rule _cluster_sv_combine_bed_filter:
    output:
        bed = CFG["dirs"]["filter"] + "bed_filter/{seq_type}--{genome_build}/combined_filter.bed"
    params:
        bed_list = CFG["options"]["filter"]["bed_list"]
    shell:
        op.as_one_line("""
        cat {params.bed_list} | 
        awk 'BEGIN{{OFS="\t"}}{{ print $1,$2,$3 }}'
        > {output.bed}
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
        > {output.bedpe}
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
        bedpe = op.switch_on_column(column = "ffpe_or_frozen", 
                                    samples = CFG["samples"], 
                                    options = ffpe_switches, 
                                    match_on = "tumour")
    output:
        bedpe = CFG["dirs"]["filter"] + "final/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.bedpe"
    wildcard_constraints:
        seq_type = "genome"
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

#######################################################################
# RUN CLUSTER SV
#######################################################################
rule _cluster_sv_run:
    input:
        bedpe = str(rules._cluster_sv_symlink_filter.output.bedpe),
        cluster_sv = str(rules._cluster_sv_install.output.cluster_sv)
    output:
        clusters = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_clusters_and_footprints.tsv",
        pval = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_distance_pvals"
    log:
        stdout = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stdout.log",
        stderr = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stderr.log"
    params:
        prefix = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}",
        chr_sizes = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/references/{genome_build}.chrom_sizes", # another script from the repo
        coords = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/references/{genome_build}_centromere_and_telomere_coords.txt" # convert output to bed format
    conda:
        CFG["conda_envs"]["cluster_sv"]
    threads:
        CFG["threads"]["cluster_sv"]
    resources:
        mem_mb = CFG["mem_mb"]["cluster_sv"]
    shell:
        op.as_one_line("""
        Rscript {input.cluster_sv}
        -bedpe {input.bedpe}
        -chr {params.chr_sizes}
        -cen_telo {params.coords}
        -out {params.prefix}
        -n {threads}
        > {log.stdout} 2> {log.stderr} &&
        sed -i '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tCLUSTER_ID\tNUM_SV\tFP_ID_LOW\tFP_ID_HIGH\tFP_COORDS_LOW\tFP_COORDS_HIGH\tP_VAL' 
        {output.clusters}
        """)

rule _cluster_sv_add_id:
    input:
        tsv = str(rules._cluster_sv_run.output.clusters)
    output:
        tsv = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_clusters_and_footprints_with_id.tsv"
    shell:
        op.as_one_line("""
        sed '/^#/d' {input.tsv} |
        sed "s/$/\t{wildcards.tumour_id}\t{wildcards.normal_id}\t{wildcards.pair_status}/" | \
        sed '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tCLUSTER_ID\tNUM_SV\tFP_ID_LOW\tFP_ID_HIGH\tFP_COORDS_LOW\tFP_COORDS_HIGH\tP_VAL\tTUMOUR_ID\tNORMAL_ID\tPAIR_STATUS' \
        > {output.tsv}
        """)

def _cluster_sv_get_clusters(wildcards):
    CFG = config["lcr-modules"]["cluster_sv"]
    SAMPLES = config["lcr-modules"]["cluster_sv"]["samples"]
    RUNS = config["lcr-modules"]["cluster_sv"]["runs"]

    if wildcards.cohort in SAMPLES.cohort:
        TUMOURS = SAMPLES.loc[SAMPLES.cohort == wildcards.cohort].sample_id
        RUNS = RUNS.loc[RUNS.tumour_sample_id.isin(TUMOURS)]

    tsv_files = expand(
        [
            str(rules._cluster_sv_add_id.output.tsv)
        ],
        zip,
        seq_type=RUNS["tumour_seq_type"],
        genome_build=RUNS["tumour_genome_build"],
        tumour_id=RUNS["tumour_sample_id"],
        normal_id=RUNS["normal_sample_id"],
        pair_status=RUNS["pair_status"])

    return { 'tsv': tsv_files }


rule _cluster_sv_combine_cluster_tsv:
    input:
        unpack(_cluster_sv_get_clusters)
    output:
        tsv = CFG["dirs"]["combine"] + "{seq_type}--{genome_build}/{cohort}.sv_clusters_and_footprints_with_id.tsv"
    shell:
        op.as_one_line("""
        awk '(NR == 1) || (FNR > 1)' {input.tsv}
        > {output.tsv}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cluster_sv_output_tsv:
    input:
        clusters = str(rules._cluster_sv_add_id.output.tsv), # str(rules._cluster_sv_add_header_tsv.output.clusters),
        pval = str(rules._cluster_sv_run.output.pval)
    output:
        clusters = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_clusters_and_footprints_id.tsv",
        pval = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_distance_pvals.tsv"
    run:
        op.relative_symlink(input.clusters, output.clusters)
        op.relative_symlink(input.pval, output.pval)


# Call clusterSV pipeline
rule _cluster_sv_all:
    input:
        expand(
            [
                str(rules._cluster_sv_output_tsv.output.clusters),
                str(rules._cluster_sv_output_tsv.output.pval),
                str(rules._cluster_sv_combine_cluster_tsv.output.tsv)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            cohort=CFG["options"]["cohort"])


##### CLEANUP #####

# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
