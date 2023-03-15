#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd

# Needed for getting snapshot paths
import os
# Needed for copying contents of individual variant batch scripts to a merged sample batch_script
import shutil

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
# `CFG` is a shortcut to `config["lcr-modules"]["igv"]`
CFG = op.setup_module(
    name = "igv",
    version = "1.0",
    subdirectories = ["inputs", "batch_scripts", "igv", "snapshots", "outputs"],
)

# Rename genome builds in metadata to match up with MAFs?
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["genome_map"]["grch37"]), "grch37", inplace=True)
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["genome_map"]["hg38"]), "hg38", inplace=True)


# Define rules to be run locally when using a compute cluster
localrules:
    _igv_symlink_regions_file,
    _igv_symlink_bam,
    _igv_symlink_bai,
    _igv_symlink_maf,
    _igv_reduce_maf_cols,
    _igv_format_regions_file,
    _igv_liftover_regions,
    _igv_filter_maf,
    _igv_create_batch_script_per_variant,
    _igv_batches_to_merge,
    _igv_download_igv,
    _igv_run,
    _igv_symlink_snapshot,
    _igv_check_snapshots,


##### FUNCTIONS #####


def get_bams(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{tumour_sample_id}}.{genome_build}.bam", genome_build=metadata[(metadata.sample_id == wildcards.tumour_sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

def get_bai(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{tumour_sample_id}}.{genome_build}.bam.bai", genome_build=metadata[(metadata.sample_id == wildcards.tumour_sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

def get_maf(wildcards):
    unix_group = config["unix_group"]
    return expand(config["lcr-modules"]["igv"]["inputs"]["maf"], allow_missing=True, unix_group=unix_group)


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _igv_symlink_regions_file:
    input:
        regions_file = CFG["inputs"]["regions_file"]
    output:
        regions_file = CFG["dirs"]["inputs"] + "regions/regions_file.txt"
    run:
        op.absolute_symlink(input.regions_file, output.regions_file)

rule _igv_symlink_bam:
    input:
        bam = get_bams
    output:
        bam = CFG["dirs"]["inputs"] + "bams/{seq_type}/{tumour_sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _igv_symlink_bai:
    input:
        bai = get_bai
    output:
        bai = CFG["dirs"]["inputs"] + "bams/{seq_type}/{tumour_sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bai, output.bai)

rule _igv_symlink_maf:
    input:
        maf = get_maf
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Filter to essential columns to prevent errors in parsing with pandas
rule _igv_reduce_maf_cols:
    input:
        maf = str(rules._igv_symlink_maf.output.maf)
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.maf.temp")
    shell:
        op.as_one_line("""
        cut -f 1,5,6,7,9,10,11,13,16 {input.maf} > {output.maf}
        """)

# Prepare regions file for liftover
rule _igv_format_regions_file:
    input:
        regions = str(rules._igv_symlink_regions_file.output.regions_file)
    output:
        regions = CFG["dirs"]["inputs"] + "regions/regions_file_formatted.txt"
    params:
        regions_format = CFG["inputs"]["regions_format"],
        oncodriveclustl_params = CFG["filter_maf"]["oncodriveclustl_options"],
        regions_build = CFG["inputs"]["regions_build"]
    log:
        stdout = CFG["logs"]["inputs"] + "format_regions.stdout.log",
        stderr = CFG["logs"]["inputs"] + "format_regions.stderr.log"
    script:
        config["lcr-modules"]["igv"]["scripts"]["format_regions"]

REGIONS_FORMAT = {
    "maf": "maf",
    "oncodriveclustl": "bed",
    "hotmaps": "bed",
    "mutation_id": "bed"
}

rule _igv_liftover_regions:
    input:
        regions = str(rules._igv_format_regions_file.output.regions),
        liftover_script = CFG["scripts"]["region_liftover_script"]
    output:
        regions = CFG["dirs"]["inputs"] + "regions/regions_file_{genome_build}.crossmap.txt"
    params:
        chain_file = reference_files(CFG["liftover_regions"]["reference_chain_file"][(CFG["inputs"]["regions_build"]).replace("hg19","grch37").replace("grch38","hg38")]),
        target_reference = lambda w: config["lcr-modules"]["igv"]["liftover_regions"]["target_reference"][w.genome_build],
        regions_type = REGIONS_FORMAT[CFG["inputs"]["regions_format"].lower()],
        regions_build = CFG["inputs"]["regions_build"].replace("grch37","GRCh37").replace("hg38","GRCh38"),
        target_build = lambda w: w.genome_build.replace("grch37","GRCh37").replace("hg38", "GRCh38")
    conda:
        CFG["conda_envs"]["liftover_regions"]
    log:
        stdout = CFG["logs"]["inputs"] + "liftover_regions_{genome_build}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "liftover_regions_{genome_build}.stderr.log"
    shell:
        op.as_one_line("""
        {input.liftover_script} {input.regions} 
        {params.regions_type} {params.regions_build} {params.target_build} 
        {output.regions} {params.chain_file} 
        {params.target_reference} > {log.stdout} 2> {log.stderr}
        """)

# Pass metadata as a pandas dataframe directly from the samples value specified in config
rule _igv_filter_maf:
    input:
        maf = str(rules._igv_reduce_maf_cols.output.maf),
        regions = str(rules._igv_liftover_regions.output.regions)
    output:
        maf = CFG["dirs"]["inputs"] + "maf/filtered_maf/{seq_type}--{genome_build}/{tumour_sample_id}--{normal_sample_id}--{pair_status}.maf"
    params:
        regions_format = REGIONS_FORMAT[CFG["inputs"]["regions_format"].lower()],
        oncodriveclustl_params = CFG["filter_maf"]["oncodriveclustl_options"],
        n_snapshots = CFG["filter_maf"]["n_snapshots"] if CFG["filter_maf"]["n_snapshots"] is not None else 1000000,
        metadata = CFG["runs"]
    wildcard_constraints:
        seq_type = "[a-zA-Z]+"
    script:
        config["lcr-modules"]["igv"]["scripts"]["filter_script"]

# Create multiple batch scripts per sample for each variant
checkpoint _igv_create_batch_script_per_variant:
    input:
        filter_maf = expand(str(rules._igv_filter_maf.output.maf), zip, normal_sample_id=CFG["runs"]["normal_sample_id"], pair_status=CFG["runs"]["pair_status"], allow_missing=True)[0],
        bam_file = str(rules._igv_symlink_bam.output.bam),
        bai_file = str(rules._igv_symlink_bai.output.bai),
        regions_lifted = str(rules._igv_liftover_regions.output.regions),
        regions_formatted = str(rules._igv_format_regions_file.output.regions)
    output:
        variant_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{tumour_sample_id}.batch"
    params:
        batch_dir = config["lcr-modules"]["igv"]["dirs"]["batch_scripts"],
        snapshot_dir = config["lcr-modules"]["igv"]["dirs"]["snapshots"],
        genome_build = lambda w: w.genome_build,
        seq_type = lambda w: w.seq_type,
        padding = config["lcr-modules"]["igv"]["generate_batch_script"]["padding"],
        igv_options = config["lcr-modules"]["igv"]["generate_batch_script"]["igv_options"],
        max_height = config["lcr-modules"]["igv"]["generate_batch_script"]["max_height"]
    script:
        config["lcr-modules"]["igv"]["scripts"]["batch_script_per_variant"]

rule _igv_batches_to_merge:
    input:
        batch_script = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{chromosome}:{start_position}--{padding}--{gene}--{tumour_sample_id}.batch"
    output:
        dispatched_batch_script = CFG["dirs"]["batch_scripts"] + "dispatched_batch_scripts/{seq_type}--{genome_build}/{chromosome}:{start_position}--{padding}--{gene}--{tumour_sample_id}.batched"
    params:
        batch_script_file = str(rules._igv_create_batch_script_per_variant.output.variant_batch)
    threads: (workflow.cores / 10)
    run:
        batch_script_path = os.path.abspath(input.batch_script)
        output_file = os.path.abspath(params.batch_script_file)

        batch_script = open(batch_script_path, "r")

        with open(output_file, "a") as handle:
            for line in batch_script:
                handle.write(line)
        batch_script.close()

        output_touch = open(output.dispatched_batch_script, "w")
        output_touch.close()

def _evaluate_batches(wildcards):
    CFG = config["lcr-modules"]["igv"]
    checkpoint_output = checkpoints._igv_create_batch_script_per_variant.get(**wildcards).output.variant_batch
    maf = expand(rules._igv_filter_maf.output.maf, zip, seq_type=wildcards.seq_type, genome_build=wildcards.genome_build, tumour_sample_id=wildcards.tumour_sample_id, normal_sample_id=CFG["runs"][(CFG["runs"]["tumour_sample_id"]==wildcards.tumour_sample_id) & (CFG["runs"]["tumour_seq_type"]==wildcards.seq_type) & (CFG["runs"]["tumour_genome_build"]==wildcards.genome_build)]["normal_sample_id"], pair_status=CFG["runs"][(CFG["runs"]["tumour_sample_id"]==wildcards.tumour_sample_id) & (CFG["runs"]["tumour_seq_type"]==wildcards.seq_type) & (CFG["runs"]["tumour_genome_build"]==wildcards.genome_build)]["pair_status"])
    
    if os.path.exists(maf[0]):
        maf_table = pd.read_table(maf[0], comment="#", sep="\t")

        return expand(
            expand(
                str(rules._igv_batches_to_merge.output.dispatched_batch_script),
                zip,
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                tumour_sample_id = maf_table["Tumor_Sample_Barcode"],
                seq_type = maf_table["seq_type"],
                genome_build = maf_table["genome_build"],
                allow_missing=True
            ),
            padding=str(CFG["generate_batch_script"]["padding"])
        )
    else:
        return []

rule _igv_download_igv:
    output:
        igv_zip = CFG["dirs"]["igv"] + "IGV_2.7.2.zip",
        igv_installed = CFG["dirs"]["igv"] + "igv_2.7.2.installed"
    conda:
        CFG["conda_envs"]["wget"]
    log:
        stdout = CFG["logs"]["igv"] + "download_igv.stdout.log",
        stderr = CFG["logs"]["igv"] + "download_igv.stderr.log"
    shell:
        op.as_one_line("""
        wget -O {output.igv_zip} https://data.broadinstitute.org/igv/projects/downloads/2.7/IGV_Linux_2.7.2.zip &&
        unzip {output.igv_zip} -d $(dirname {output.igv_zip}) > {log.stdout} 2> {log.stderr} &&
        touch {output.igv_installed}
        """)

checkpoint _igv_run:
    input:
        igv = str(rules._igv_download_igv.output.igv_installed),
        batch_script = _evaluate_batches
    output:
        complete = CFG["dirs"]["snapshots"] + "completed/{seq_type}--{genome_build}--{tumour_sample_id}.completed"
    log:
        stdout = CFG["logs"]["igv"] + "{seq_type}--{genome_build}/{tumour_sample_id}_igv_run.stdout",
        stderr = CFG["logs"]["igv"] + "{seq_type}--{genome_build}/{tumour_sample_id}_igv_run.stderr"
    params:
        merged_batch = str(rules._igv_create_batch_script_per_variant.output.variant_batch),
        igv = CFG["dirs"]["igv"] + "IGV_Linux_2.7.2/igv.sh"
    threads: (workflow.cores / 5)
    shell:
        op.as_one_line("""
        echo 'exit' >> {params.merged_batch} ;
        xvfb-run --auto-servernum {params.igv} -b {params.merged_batch} > {log.stdout} 2> {log.stderr} ;
        touch {output.complete}
        """)

rule _igv_symlink_snapshot:
    input:
        snapshot = CFG["dirs"]["snapshots"] + "{seq_type}--{genome_build}/{chromosome}/{chromosome}:{start_position}--{padding}--{gene}--{tumour_sample_id}.png"
    output:
        snapshot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{chromosome}/{chromosome}:{start_position}--{padding}--{gene}--{tumour_sample_id}.png"
    run:
        op.relative_symlink(input.snapshot, output.snapshot)

def _symlink_snapshot(wildcards):
    CFG = config["lcr-modules"]["igv"]
    checkpoint_outputs = checkpoints._igv_run.get(**wildcards).output.complete
    maf = expand(
        str(rules._igv_filter_maf.output.maf), 
        zip, seq_type=wildcards.seq_type, 
        genome_build=wildcards.genome_build, 
        tumour_sample_id=wildcards.tumour_sample_id, 
        normal_sample_id=CFG["runs"][(CFG["runs"]["tumour_sample_id"]==wildcards.tumour_sample_id) & (CFG["runs"]["tumour_seq_type"]==wildcards.seq_type) & (CFG["runs"]["tumour_genome_build"]==wildcards.genome_build)]["normal_sample_id"], 
        pair_status=CFG["runs"][(CFG["runs"]["tumour_sample_id"]==wildcards.tumour_sample_id) & (CFG["runs"]["tumour_seq_type"]==wildcards.seq_type) & (CFG["runs"]["tumour_genome_build"]==wildcards.genome_build)]["pair_status"]
    )

    if os.path.exists(maf[0]):
        maf_table = pd.read_table(maf[0], comment="#", sep="\t")

        return expand(
            expand(
                str(rules._igv_symlink_snapshot.output.snapshot), 
                zip,
                seq_type = maf_table["seq_type"],
                genome_build = maf_table["genome_build"],
                chromosome = maf_table["chr_std"],
                start_position = maf_table["Start_Position"],
                gene = maf_table["Hugo_Symbol"],
                tumour_sample_id = maf_table["Tumor_Sample_Barcode"],
                allow_missing=True
            ),
            padding = str(CFG["generate_batch_script"]["padding"])
        )
    else:
        return []

rule _igv_check_snapshots:
    input:
        snapshots = _symlink_snapshot,
        igv_completed = str(rules._igv_run.output.complete)
    output:
        snapshots = CFG["dirs"]["outputs"] + "completed/{seq_type}--{genome_build}--{tumour_sample_id}.completed"
    shell:
        "touch {output.snapshots}"

# Generates the target sentinels for each run, which generate the symlinks
if CFG["test_run"] is False:
    rule _igv_all:
        input:
            expand([str(rules._igv_run.output.complete), str(rules._igv_check_snapshots.output.snapshots)], zip, tumour_sample_id=CFG["runs"]["tumour_sample_id"], seq_type=CFG["runs"]["tumour_seq_type"], genome_build=CFG["runs"]["tumour_genome_build"])

#if CFG["test_run"] is False:
#    rule _igv_all:
#        input:
#            expand(str(rules._igv_check_outputs.output), zip, tumour_sample_id=CFG["runs"]["tumour_sample_id"], seq_type=CFG["runs"]["tumour_seq_type"], genome_build=CFG["runs"]["tumour_genome_build"])

if CFG["test_run"] is True:
    rule _igv_all:
        input:
            expand(rules._igv_filter_maf.output.maf, zip, seq_type=CFG["runs"]["tumour_seq_type"], tumour_sample_id=CFG["runs"]["tumour_sample_id"], normal_sample_id=CFG["runs"]["normal_sample_id"], pair_status=CFG["runs"]["pair_status"], genome_build=CFG["runs"]["tumour_genome_build"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
