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
    subdirectories = ["inputs", "reformat_bedpe", "cluster_sv", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_install,
    _cluster_sv_input_bedpe,
    _cluster_sv_reformat_bedpe,
    _cluster_sv_output_tsv,
    _cluster_sv_all,


##### RULES #####

# download cluster files from git repo without cloning the repo itself
# decompress files into the 00-inputs

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cluster_sv_input_bedpe:
    input:
        bedpe = CFG["inputs"]["sample_bedpe"]
    output:
        bedpe = CFG["dirs"]["inputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)

CFG = op.setup_module(
    name = "cluster_sv",
    version = "1.0",
    subdirectories = ["inputs", "filter_vcf", "annotate", "FFPE_filter",
                      "combine", "clusterSV", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_install,
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
        vcf = temp(CFG["dirs"]["filter_vcf"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf")
    shell:
        "gzip -dc {input.vcf_gz} > {output.vcf}" #this should work on both gzip and bcftools compressed files 


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _cluster_sv_filter_vcf:
    input:
        vcf = str(rules._cluster_sv_decompress_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["filter_vcf"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.filtered.vcf"
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


def _cluster_sv_get_FFPE_bedpe(wildcards):
    CFG = config["lcr-modules"]["cluster_sv"]
    RUNS = config["lcr-modules"]["cluster_sv"]["runs"]

    if CFG["options"]["FFPE_filter"] is None:
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

rule _cluster_sv_FFPE_small_INV_filter:
    input:
        bedpe = str(rules._cluster_sv_annotate.output.bedpe)
    output:
        bedpe = CFG["dirs"]["filtered"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.annotated.bedpe"
    log: 
        CFG["logs"]["annotate"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_R_annotateSV.log"
    params:
        r_script = CFG["options"]["vcf_annotation_script"]
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


rule _cluster_sv_run:
    input:
        bedpe = str(rules._cluster_sv_reformat_bedpe.output.bedpe),
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


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cluster_sv_output_tsv:
    input:
        clusters = str(rules._cluster_sv_run.output.clusters),#str(rules._cluster_sv_add_header_tsv.output.clusters),
        pval = str(rules._cluster_sv_run.output.pval)
    output:
        clusters = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_clusters_and_footprints.tsv",
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
                str(rules._cluster_sv_output_tsv.output.pval)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


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

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
