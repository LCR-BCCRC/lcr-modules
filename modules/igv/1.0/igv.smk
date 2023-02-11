#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd

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
CFG["samples"]["genome_build"].mask(CFG["samples"]["genome_build"].isin(CFG["genome_map"]["grch37"]), "grch37", inplace=True)
CFG["samples"]["genome_build"].mask(CFG["samples"]["genome_build"].isin(CFG["genome_map"]["hg38"]), "hg38", inplace=True)


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _igv_symlink_regions_file,
    _igv_symlink_metadata,
    _igv_symlink_maf,
    _igv_format_regions_file,
    _igv_liftover_regions,
    _igv_create_batch_script,
    _igv_download_igv,
    _igv_run

##### FUNCTIONS #####


def get_bams(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{sample_id}}.{genome_build}.bam", genome_build=metadata[(metadata.sample_id == wildcards.sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])



##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _igv_symlink_regions_file:
    input:
        regions_file = CFG["inputs"]["regions_file"]
    output:
        regions_file = CFG["dirs"]["inputs"] + "regions/regions_file.txt"
    run:
        op.absolute_symlink(input.regions_file, output.regions_file)

rule _igv_symlink_bams:
    input:
        bam = get_bams
    output:
        bam = CFG["dirs"]["inputs"] + "bams/{seq_type}/{sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _igv_symlink_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

rule _igv_reduce_maf_cols:
    input:
        maf = str(rules._igv_symlink_maf.output.maf)
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}_cols.maf")
    shell:
        op.as_one_line("""
        cut -f 1,5,6,7,9,10,11,13,16 {input.maf} > {output.maf}
        """)

# Prepare regions file for liftover
rule _igv_format_regions_file:
    input:
        regions = str(rules._igv_symlink_regions_file.output.regions_file)
    output:
        regions = config["lcr-modules"]["igv"]["dirs"]["inputs"] + "regions/regions_file_formatted.txt"
    params:
        regions_format = config["lcr-modules"]["igv"]["inputs"]["regions_format"],
        oncodriveclustl_params = config["lcr-modules"]["igv"]["filter_maf"]["oncodriveclustl_options"],
        regions_build = config["lcr-modules"]["igv"]["inputs"]["regions_build"]
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
        regions_lifted = CFG["dirs"]["inputs"] + "regions/regions_file_{genome_build}.txt"
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
        {output.regions_lifted} {params.chain_file} 
        {params.target_reference} > {log.stdout} 2> {log.stderr}
        """)

# Pass metadata as a pandas dataframe directly from the samples value specified in config
rule _igv_filter_maf:
    input:
        maf = str(rules._igv_reduce_maf_cols.output.maf),
        regions = str(rules._igv_liftover_regions.output.regions_lifted)
    output:
        maf_filtered = config["lcr-modules"]["igv"]["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}_cols_filtered.maf"
    params:
        regions_format = REGIONS_FORMAT[config["lcr-modules"]["igv"]["inputs"]["regions_format"].lower()],
        metadata = config["lcr-modules"]["igv"]["samples"],
        genome_build = lambda w: w.genome_build,
        seq_type = lambda w: w.seq_type,
        genome_map = config["lcr-modules"]["igv"]["genome_map"],
        oncodriveclustl_params = config["lcr-modules"]["igv"]["filter_maf"]["oncodriveclustl_options"],
        n_snapshots = config["lcr-modules"]["igv"]["filter_maf"]["n_snapshots"] if config["lcr-modules"]["igv"]["filter_maf"]["n_snapshots"] is not None else 1000000
    script:
        config["lcr-modules"]["igv"]["scripts"]["filter_script"]
    
# Pass filtered MAF to create batch script
rule _igv_create_batch_script:
    input:
        maf_filtered = str(rules._igv_filter_maf.output.maf_filtered)
    output:
        batch_script = temp(config["lcr-modules"]["igv"]["dirs"]["batch_scripts"] + "{seq_type}--{genome_build}.batch")
    params:
        metadata = config["lcr-modules"]["igv"]["samples"],
        snapshot_dir = config["lcr-modules"]["igv"]["dirs"]["snapshots"],
        genome_build = lambda w: w.genome_build,
        seq_type = lambda w: w.seq_type,
        padding = config["lcr-modules"]["igv"]["generate_batch_script"]["padding"],
        max_height = config["lcr-modules"]["igv"]["generate_batch_script"]["max_height"]
    script:
        config["lcr-modules"]["igv"]["scripts"]["batch_script"]

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

rule _igv_track_samples:
    input:
        bam = str(rules._igv_symlink_bams.output.bam)
    output:
        finished = CFG["dirs"]["outputs"] + "samples/{seq_type}/{sample_id}_" + CFG["inputs"]["regions_format"] + ".track"
    shell:
        "touch {output.finished}"

rule _igv_run:
    input:
        batch_script = str(rules._igv_create_batch_script.output.batch_script),
        igv_installed = str(rules._igv_download_igv.output.igv_installed),
        sample_track = expand(str(rules._igv_track_samples.output.finished), zip, seq_type=CFG["samples"]["seq_type"], sample_id=CFG["samples"]["sample_id"])
    output:
        success = CFG["dirs"]["outputs"] + "snapshots/{seq_type}--{genome_build}_snapshot.finished"
    params:
        #igv = CFG["dirs"]["igv"] + "IGV_Linux_2.7.2/igv.sh"
        igv = "/projects/rmorin/projects/RNA_seq_ssm/test/bin/IGV_Linux_2.7.2/igv.sh"
    log:
        stdout = CFG["logs"]["igv"] + "run_igv_{seq_type}--{genome_build}.stdout.log",
        stderr = CFG["logs"]["igv"] + "run_igv_{seq_type}--{genome_build}.stderr.log"
    shell:
        op.as_one_line("""
        xvfb-run --auto-servernum {params.igv} -b {input.batch_script} > {log.stdout} 2> {log.stderr} &&
        touch {output.success}
        """)

# Generates the target sentinels for each run, which generate the symlinks
if CFG["test_run"] is False:
    rule _igv_all:
        input:
            expand(str(rules._igv_run.output.success), seq_type=CFG["samples"]["seq_type"], genome_build=CFG["samples"]["genome_build"])

if CFG["test_run"] is True:
    rule _igv_all:
        input:
            expand(rules._igv_filter_maf.output.maf_filtered, seq_type=CFG["samples"]["seq_type"], genome_build=CFG["samples"]["seq_type"])




##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
