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


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _igv_symlink_regions_file,
    _igv_symlink_metadata,
    _igv_symlink_maf,
    _igv_liftover_regions,
    _igv_create_batch_script,
    _igv_download_igv,
    _igv_run


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _igv_symlink_regions_file:
    input:
        regions_file = CFG["inputs"]["regions_file"]
    output:
        regions_file = CFG["dirs"]["inputs"] + "regions/regions_file.txt"
    run:
        op.absolute_symlink(input.regions_file, output.regions_file)

rule _igv_symlink_metadata:
    input:
        metadata = CFG["inputs"]["metadata"]
    output:
        metadata = CFG["dirs"]["inputs"] + "metadata/metadata.tsv"
    run:
        op.absolute_symlink(input.metadata, output.metadata)

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


rule _igv_liftover_regions:
    input:
        regions = str(rules._igv_symlink_regions_file.output.regions_file),
        liftover_script = CFG["scripts"]["region_liftover_script"]
    output:
        regions_lifted = CFG["dirs"]["inputs"] + "regions/regions_file_{genome_build}.txt"
    params:
        chain_file = reference_files(CFG["liftover_regions"]["reference_chain_file"][(CFG["inputs"]["regions_build"]).replace("hg19","grch37").replace("grch38","hg38")]),
        target_reference = lambda w: config["lcr-modules"]["igv"]["liftover_regions"]["target_reference"][w.genome_build],
        regions_type = CFG["inputs"]["regions_format"].lower(),
        target_build = lambda w: w.genome_build.replace("grch37","GRCh37").replace("hg38", "GRCh38")
    conda:
        CFG["conda_envs"]["liftover_regions"]
    log:
        stdout = CFG["logs"]["inputs"] + "liftover_regions_{genome_build}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "liftover_regions_{genome_build}.stderr.log"
    shell:
        op.as_one_line("""
        {input.liftover_script} {input.regions} 
        {params.regions_type} {params.target_build} 
        {output.regions_lifted} {params.chain_file} 
        {params.target_reference} > {log.stdout} 2> {log.stderr}
        """)

# Filter the MAF based on regions file
rule _igv_filter_maf:
    input:
        maf = str(rules._igv_reduce_maf_cols.output.maf),
        regions = str(rules._igv_liftover_regions.output.regions_lifted)
    output:
        maf_filtered = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}_cols_filtered.maf"
    params:
        regions_format = CFG["inputs"]["regions_format"].lower(),
        metadata = CFG["inputs"]["metadata"],
        genome_build = lambda w: w.genome_build,
        seq_type = lambda w: w.seq_type,
        genome_map = CFG["genome_map"],
        oncodriveclustl_params = CFG["filter_maf"]["oncodriveclustl_options"]
    run:
        # Read input MAF and regions into pandas df
        maf_df = pd.read_table(input.maf, comment="#", sep="\t")
        regions_df = pd.read_table(input.regions, comment="#", sep="\t")
        if params.regions_format in ["maf"]:
        # Create common columns to subset the larger MAF down
            count = 0
            for df in [maf_df, regions_df]:
                count += 1
                if count == 1:
                    print(f"Working on maf df")
                if count == 2:
                    print(f"Working on regions df")
                df["chr_std"] = df.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
                df["genomic_pos_std"] = df["chr_std"] + ":" + df["Start_Position"].map(str) + "_" + df["End_Position"].map(str)

            # Filter larger MAF 
            filtered_maf = maf_df[maf_df["genomic_pos_std"].isin(regions_df["genomic_pos_std"])]

            # Filter only to BAM files of corresponding build and seq_type
            SAMPLES = op.load_samples(params.metadata)
            SAMPLES = op.filter_samples(SAMPLES, seq_type=params.seq_type)
            genome_build_list = params.genome_map[params.genome_build]
            print(f"Only including samples of these builds: {genome_build_list}")
            BUILD_SAMPLES = op.filter_samples(SAMPLES, genome_build=genome_build_list)
            filtered_maf = filtered_maf[filtered_maf["Tumor_Sample_Barcode"].isin(BUILD_SAMPLES.sample_id)]

            # Write output
            filtered_maf.to_csv(output.maf_filtered, sep="\t")

    
# Pass filtered MAF to create batch script
rule _igv_create_batch_script:
    input:
        maf_filtered = str(rules._igv_filter_maf.output.maf_filtered),
        metadata = str(rules._igv_symlink_metadata.output.metadata)
    output:
        batch_script = temp(CFG["dirs"]["batch_scripts"] + "{seq_type}--{genome_build}.batch")
    params:
        py_script = CFG["scripts"]["batch_script"],
        snapshot_dir = CFG["dirs"]["snapshots"],
        genome_build = lambda w: w.genome_build,
        seq_type = lambda w: w.seq_type,
        padding = CFG["generate_batch_script"]["padding"],
        max_height = CFG["generate_batch_script"]["max_height"],
        n_snapshots = CFG["generate_batch_script"]["n_snapshots"]
    log:
        stdout = CFG["logs"]["batch_scripts"] + "{seq_type}--{genome_build}_batch_script.stdout.log",
        stderr = CFG["logs"]["batch_scripts"] + "{seq_type}--{genome_build}_batch_script.stderr.log"
    shell:
        op.as_one_line("""
        {params.py_script} {input.maf_filtered} 
        --output {output.batch_script} --metadata {input.metadata} 
        --padding {params.padding} --max_height {params.max_height} 
        --snapshot_dir {params.snapshot_dir} --n_snapshots {params.n_snapshots}
        --genome_build {params.genome_build} --seq_type {params.seq_type} > {log.stdout} 2> {log.stderr}
        """)

#rule _igv_merge_batch_scripts:
#    input:
#        batch_scripts = expand(str(rules._igv_create_batch_script.output.batch_script), genome_build=["hg38","grch37"], seq_type=["capture","genome"]),
#    output:
#        merged_batch = CFG["dirs"]["batch_scripts"] + "merged_script.batch"
#    params:
#        script_dir = CFG["dirs"]["batch_scripts"]
#    shell:
#        op.as_one_line("""
#        batch_dir={params.script_dir} &&
#        cat <(cat $(echo $batch_dir)/*.batch) <(echo end) | awk '{{ if ($0 !~ "exit") print $0 }}' | sed 's/end/exit\n/g' > {output.merged_batch}
#        """)


#### WHEN LAST CHECKED DRY RUN WORKS UP TO HERE B-) 
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

rule _igv_run:
    input:
        batch_script = str(rules._igv_create_batch_script.output.batch_script),
        igv_installed = str(rules._igv_download_igv.output.igv_installed),
    output:
        success = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}_snapshot.finished"
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

#rule _igv_run:
#    input:
#        batch_script = str(rules._igv_merge_batch_scripts.output.merged_batch),
#        igv_installed = str(rules._igv_download_igv.output.igv_installed),
#    output:
#        success = CFG["dirs"]["outputs"] + "merged_batch_snapshot.finished"
#    params:
#        igv = CFG["dirs"]["igv"] + "IGV_Linux_2.16.0/igv.sh"
#    log:
#        stdout = CFG["logs"]["igv"] + "run_igv_merged_batch.stdout.log",
#        stderr = CFG["logs"]["igv"] + "run_igv_merged_batch.stderr.log"
#    shell:
#        op.as_one_line("""
#        xvfb-run --auto-servernum {params.igv} -b {input.batch_script} > {log.stdout} 2> {log.stderr} &&
#        touch {output.success}
#        """)

# Generates the target sentinels for each run, which generate the symlinks
rule _igv_all:
    input:
        expand(rules._igv_run.output.success, seq_type=["genome","capture"], genome_build=["hg38","grch37"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
