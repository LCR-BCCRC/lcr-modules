#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd
import os

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

# Rename genome_build values in sample metadata to correlate with MAF values
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["options"]["genome_map"]["grch37"]), "grch37", inplace=True)
CFG["runs"]["tumour_genome_build"].mask(CFG["runs"]["tumour_genome_build"].isin(CFG["options"]["genome_map"]["hg38"]), "hg38", inplace=True)

# Define output file suffix based on config parameters
SUFFIX = ".pad" + str(CFG["options"]["generate_batch_script"]["padding"])

# Assign pair_status_directory value based on pair_status value
PAIR_STATUS_DICT = {"matched": "tumour_normal_pair", "unmatched": "tumour_only"}

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
    _igv_touch_summary,
    _igv_estimate_snapshots,
    _igv_snapshot_estimate_finished


##### FUNCTIONS #####


def get_bams(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{tumour_id}}.{genome_build}.bam", genome_build=metadata[(metadata.sample_id == wildcards.tumour_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

def get_bai(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{tumour_id}}.{genome_build}.bam.bai", genome_build=metadata[(metadata.sample_id == wildcards.tumour_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

def get_normal_bam(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{normal_sample_id}}.{genome_build}.bam", genome_build=metadata[(metadata.sample_id == wildcards.normal_sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

def get_normal_bai(wildcards):
    metadata = config["lcr-modules"]["igv"]["samples"]
    return expand("data/{{seq_type}}_bams/{{normal_sample_id}}.{genome_build}.bam.bai", genome_build=metadata[(metadata.sample_id == wildcards.normal_sample_id) & (metadata.seq_type == wildcards.seq_type)]["genome_build"])

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
        bam = CFG["dirs"]["inputs"] + "bams/{seq_type}/{tumour_id}.bam"
    threads:
        CFG["threads"]["_igv_symlink_bam"]
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _igv_symlink_bai:
    input:
        bai = get_bai
    output:
        bai = CFG["dirs"]["inputs"] + "bams/{seq_type}/{tumour_id}.bam.bai"
    threads:
        CFG["threads"]["_igv_symlink_bai"]
    run:
        op.absolute_symlink(input.bai, output.bai)

rule _igv_symlink_normal_bam:
    input:
        bam = get_normal_bam
    output:
        bam = CFG["dirs"]["inputs"] + "normal_bams/{seq_type}/{normal_sample_id}.bam"
    run:
        op.absolute_symlink(input.bam, output.bam)

rule _igv_symlink_normal_bai:
    input:
        bai = get_normal_bai
    output:
        bai = CFG["dirs"]["inputs"] + "normal_bams/{seq_type}/{normal_sample_id}.bam.bai"
    run:
        op.absolute_symlink(input.bai, output.bai)

rule _igv_symlink_maf:
    input:
        maf = get_maf
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Reduce MAF columns to prevent parsing errors in Pandas
rule _igv_reduce_maf_cols:
    input:
        maf = str(rules._igv_symlink_maf.output.maf)
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf.temp")
    shell:
        op.as_one_line("""
        cut -f 1,5,6,7,9,10,11,13,16 {input.maf} > {output.maf}
        """)

# Convert input regions file into BED format
rule _igv_format_regions_file:
    input:
        regions = str(rules._igv_symlink_regions_file.output.regions_file)
    output:
        regions = CFG["dirs"]["inputs"] + "regions/regions_file_formatted.txt"
    params:
        regions_format = CFG["inputs"]["regions_format"],
        oncodriveclustl_params = CFG["options"]["filter_maf"]["oncodriveclustl_options"],
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
        chain_file = reference_files(CFG["options"]["liftover_regions"]["reference_chain_file"][(CFG["inputs"]["regions_build"]).replace("hg19","grch37").replace("grch38","hg38")]),
        target_reference = lambda w: config["lcr-modules"]["igv"]["options"]["liftover_regions"]["target_reference"][w.genome_build],
        regions_type = REGIONS_FORMAT[CFG["inputs"]["regions_format"].lower()],
        regions_build = CFG["inputs"]["regions_build"].replace("grch37","GRCh37").replace("hg38","GRCh38"),
        target_build = lambda w: w.genome_build.replace("grch37","GRCh37").replace("hg38", "GRCh38")
    conda:
        CFG["conda_envs"]["liftover_regions"]
    resources:
        **CFG["resources"]["_igv_liftover_regions"]
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

# Filter MAF to lines containing positions of interest
rule _igv_filter_maf:
    input:
        maf = str(rules._igv_reduce_maf_cols.output.maf),
        regions = str(rules._igv_liftover_regions.output.regions)
    output:
        maf = CFG["dirs"]["inputs"] + "maf/filtered_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}.maf"
    params:
        regions_format = REGIONS_FORMAT[CFG["inputs"]["regions_format"].lower()],
        oncodriveclustl_params = CFG["options"]["filter_maf"]["oncodriveclustl_options"],
        metadata = CFG["runs"]
    log:
        stdout = CFG["logs"]["inputs"] + "filter_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}/filter_maf.stdout.log",
        stderr = CFG["logs"]["inputs"] + "filter_maf/{seq_type}--{genome_build}/{tumour_id}--{normal_sample_id}--{pair_status}/filter_maf.stderr.log"
    wildcard_constraints:
        seq_type = "[a-zA-Z]+"
    script:
        config["lcr-modules"]["igv"]["scripts"]["filter_script"]

if CFG["estimate_only"] == False and CFG["identify_failed_snaps"]==False:

    def _get_maf(wildcards):
        CFG = config["lcr-modules"]["igv"]

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.tumour_id, tumour_seq_type=wildcards.seq_type)
        genome_build = this_sample["tumour_genome_build"]
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        return (
            expand(
                str(rules._igv_filter_maf.output.maf),
                zip,
                seq_type = wildcards.seq_type,
                genome_build = genome_build,
                tumour_id = wildcards.tumour_id,
                normal_sample_id = normal_sample_id,
                pair_status = pair_status
            )
        )

    def _get_bam_files(wildcards):
        CFG = config["lcr-modules"]["igv"]

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.tumour_id, tumour_seq_type=wildcards.seq_type)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        if pair_status.item() == "matched":
            tumour_bam_file = expand(str(rules._igv_symlink_bam.output.bam), zip, seq_type=wildcards.seq_type, tumour_id=wildcards.tumour_id)[0]
            normal_bam_file = expand(str(rules._igv_symlink_normal_bam.output.bam), zip, seq_type=wildcards.seq_type, normal_sample_id=normal_sample_id)[0]
            return([tumour_bam_file, normal_bam_file])

        if pair_status.item() == "unmatched":
            return (
                expand(str(rules._igv_symlink_bam.output.bam), zip, seq_type=wildcards.seq_type, tumour_id=wildcards.tumour_id)
            )

    def _get_bai_files(wildcards):
        CFG = config["lcr-modules"]["igv"]

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id=wildcards.tumour_id, tumour_seq_type=wildcards.seq_type)
        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        if pair_status.item() == "matched":
            tumour_bai_file = expand(str(rules._igv_symlink_bai.output.bai), zip, seq_type=wildcards.seq_type, tumour_id=wildcards.tumour_id)[0]
            normal_bai_file = expand(str(rules._igv_symlink_normal_bai.output.bai), zip, seq_type=wildcards.seq_type, normal_sample_id=normal_sample_id)[0]
            return([tumour_bai_file, normal_bai_file])
        if pair_status.item() == "unmatched":
            return(
                expand(str(rules._igv_symlink_bai.output.bai), zip, seq_type=wildcards.seq_type, tumour_id=wildcards.tumour_id)
            )

    # Create batch scripts for each variant
    checkpoint _igv_create_batch_script_per_variant:
        input:
            bam_file = _get_bam_files,
            bai_file = _get_bai_files,
            filter_maf = _get_maf,
            regions_lifted = str(rules._igv_liftover_regions.output.regions),
            regions_formatted = str(rules._igv_format_regions_file.output.regions)
        output:
            batches_finished = CFG["dirs"]["batch_scripts"] + "completed/{seq_type}--{genome_build}/{tumour_id}.finished",
            variant_batch = expand(CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{{seq_type}}--{{genome_build}}/{preset}/{{tumour_id}}" + SUFFIX + ".batch", preset = CFG["options"]["igv_presets"], allow_missing=True),
        params:
            batch_dir = config["lcr-modules"]["igv"]["dirs"]["batch_scripts"],
            snapshot_dir = config["lcr-modules"]["igv"]["dirs"]["snapshots"],
            genome_build = lambda w: w.genome_build,
            seq_type = lambda w: w.seq_type,
            batch_options = config["lcr-modules"]["igv"]["options"]["generate_batch_script"],
            suffix = SUFFIX,
            igv_presets = config["lcr-modules"]["igv"]["options"]["igv_presets"]
        log:
            stdout = CFG["logs"]["batch_scripts"] + "_igv_create_batch_script_per_variant/{seq_type}--{genome_build}/{tumour_id}" + SUFFIX + ".stdout.log",
            stderr = CFG["logs"]["batch_scripts"] + "_igv_create_batch_script_per_variant/{seq_type}--{genome_build}/{tumour_id}" + SUFFIX + ".stderr.log"
        script:
            config["lcr-modules"]["igv"]["scripts"]["batch_script_per_variant"]

    # Keep track of which variant and sample_id combinations have been seen, merge individual variant batch scripts into a large batch script per sample_id
    rule _igv_batches_to_merge:
        input:
            batch_script = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".batch"
        output:
            dispatched_batch_script = CFG["dirs"]["batch_scripts"] + "dispatched_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".batch"
        params:
            merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{tumour_id}" + SUFFIX + ".batch",
            igv_options = lambda w: config["lcr-modules"]["igv"]["options"]["generate_batch_script"]["igv_options"][w.preset]
        threads: (workflow.cores / 10)
        run:
            batch_script_path = os.path.abspath(input.batch_script)
            output_file = os.path.abspath(params.merged_batch)

            batch_script = open(batch_script_path, "r")

            with open(output_file, "r") as f:
                merged_lines = len(f.readlines())

            with open(output_file, "a") as handle:
                for line in batch_script:
                    if merged_lines > 0:
                        if line.startswith(("load","maxPanelHeight","genome", "setSleepInterval", "collapse")):
                            continue
                        if line.startswith(tuple(params.igv_options)) and not line.startswith("sort"):
                            continue
                    handle.write(line)
            batch_script.close()

            output_touch = open(output.dispatched_batch_script, "w")
            output_touch.close()

    # Return list of all batch scripts that were created from the filtered maf and merged
    def _evaluate_batches(wildcards):
        CFG = config["lcr-modules"]["igv"]
        checkpoint_output = checkpoints._igv_create_batch_script_per_variant.get(**wildcards).output.variant_batch

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.tumour_id, tumour_genome_build = wildcards.genome_build, tumour_seq_type = wildcards.seq_type)

        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf), 
            zip, 
            seq_type=wildcards.seq_type, 
            genome_build=wildcards.genome_build, 
            tumour_id=wildcards.tumour_id, 
            normal_sample_id=normal_sample_id, 
            pair_status=pair_status
        )

        if os.path.exists(maf[0]):
            maf_table = pd.read_table(maf[0], comment="#", sep="\t")

            return expand(
                    expand(
                        str(rules._igv_batches_to_merge.output.dispatched_batch_script),
                        zip,
                        chromosome = maf_table["chr_std"],
                        start_position = maf_table["Start_Position"],
                        gene = maf_table["Hugo_Symbol"],
                        tumour_id = maf_table["Tumor_Sample_Barcode"],
                        seq_type = maf_table["seq_type"],
                        genome_build = maf_table["genome_build"],
                        allow_missing = True
                    ),
                    preset = wildcards.preset
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

    # Run IGV once all individual variant batch scripts have been merged into one script per sample_id
    checkpoint _igv_run:
        input:
            igv = str(rules._igv_download_igv.output.igv_installed),
            batch_script = _evaluate_batches,
            finished_batches = str(rules._igv_create_batch_script_per_variant.output.batches_finished)
        output:
            complete = CFG["dirs"]["snapshots"] + "completed/{seq_type}--{genome_build}--{tumour_id}--{preset}.completed"
        params:
            merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{tumour_id}" + SUFFIX + ".batch",
            igv = CFG["dirs"]["igv"] + "IGV_Linux_2.7.2/igv.sh",
            max_time = CFG["options"]["generate_batch_script"]["sleep_timer"],
            server_number = "-n " + CFG["options"]["xvfb_parameters"]["server_number"] if CFG["options"]["xvfb_parameters"]["server_number"] is not None else "--auto-servernum",
            server_args = CFG["options"]["xvfb_parameters"]["server_args"]
        resources:
            **CFG["resources"]["_igv_run"]
        log:
            stdout = CFG["logs"]["igv"] + "{seq_type}--{genome_build}/{preset}/{tumour_id}_igv_run.stdout.log",
            stderr = CFG["logs"]["igv"] + "{seq_type}--{genome_build}/{preset}/{tumour_id}_igv_run.stderr.log"
        threads: (workflow.cores)
        shell:
            op.as_one_line("""
            lines=$(wc -l < {params.merged_batch}) ;
            if [ $lines -gt 0 ] ;
            then
            if ! grep -q -e "exit" {params.merged_batch} ;
            then
            echo 'exit' >> {params.merged_batch} ;
            fi ;
            maxtime=$(($(wc -l < {params.merged_batch}) * 60 + 15)) ;
            timeout $maxtime xvfb-run -s "-screen 0 1920x1080x24" {params.server_number} {params.server_args} {params.igv} -b {params.merged_batch} > {log.stdout} 2> {log.stderr} && touch {output.complete} ;
            exit=$? ;
            if [ $exit -ne 0 ] ;
            then
            if grep -q -e "No such process" {log.stderr} && grep -q -e "Executing Command: exit" {log.stdout} ;
            then
            echo "All IGV batch script commands have completed succesfully, but an Xvfb-run kill error has occurred." >> {log.stdout} && touch {output.complete} ;
            else
            false ;
            fi ;
            fi ;
            else
            echo 'Skipping sample {wildcards.tumour_id} because it either has no variants to snapshot or all variants have been already been snapshot.' ;
            touch {output.complete} ;
            fi
            """)

    rule _igv_track_failed:
        output:
            failed_summary = CFG["dirs"]["outputs"] + "snapshot_estimates/failed_summary.txt",
            ready = temp(CFG["dirs"]["outputs"] + "snapshot_estimates/failed_summary.completed")
        run:
            header = "\t".join(["sample_id","seq_type","genome_build","gene","chromosome","position","preset","snapshot_path"])
            with open(output.failed_summary, "w") as handle:
                handle.write(header + "\n")
            ready = open(output.ready, "w")
            ready.close()

    rule _igv_quality_control:
        input:
            igv = str(rules._igv_run.output.complete),
            snapshot = CFG["dirs"]["snapshots"] + "{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".png",
            failed_summary = str(rules._igv_track_failed.output.ready)
        output:
            snapshot_qc = temp(CFG["dirs"]["snapshots"] + "qc/{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".qc")
        params:
            batch_script = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".batch",
            merged_batch = CFG["dirs"]["batch_scripts"] + "merged_batch_scripts/{seq_type}--{genome_build}/{preset}/{tumour_id}" + SUFFIX + ".batch",
            igv = CFG["dirs"]["igv"] + "IGV_Linux_2.7.2/igv.sh",
            server_number = "-n " + CFG["options"]["xvfb_parameters"]["server_number"] if CFG["options"]["xvfb_parameters"]["server_number"] is not None else "--auto-servernum",
            server_args = CFG["options"]["xvfb_parameters"]["server_args"],
            batch_temp = CFG["dirs"]["batch_scripts"] + "single_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".batch.temp",
            failed_summary = str(rules._igv_track_failed.output.failed_summary)
        resources:
            **CFG["resources"]["_igv_quality_control"]
        log:
            stdout = CFG["logs"]["snapshots"] + "{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + "_quality_control.stdout.log",
            stderr = CFG["logs"]["snapshots"] + "{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + "_quality_control.stderr.log"
        threads: (workflow.cores)
        run:
            import subprocess
            success = True
            corrupt_checks = 0
            is_corrupt = True
            height = None
            width = None
            while corrupt_checks < 2 and is_corrupt == True:
                corrupt_checks += 1
                try:
                    height = str(subprocess.check_output(f"identify -format '%h' {input.snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                    width = str(subprocess.check_output(f"identify -format '%w' {input.snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                    is_corrupt = False
                except:
                    os.system(f'echo "Snapshot may be corrupt. Rerunning IGV, attempt {str(corrupt_checks)}:" >> {log.stdout}')
                    os.system(f'cat {params.batch_script} > {params.batch_temp} && echo "exit" >> {params.batch_temp}')
                    os.system(f'maxtime=$(($(wc -l < {params.batch_temp}) * 60 + 15)) ; timeout --foreground $maxtime xvfb-run -s "-screen 0 1980x1020x24" {params.server_number} {params.server_args} {params.igv} -b {params.batch_temp} >> {log.stdout} 2>> {log.stderr}')
            if is_corrupt == True:
                os.system(f'echo "Snapshot may be corrupt." >> {log.stdout}')
                success = False
            if height in ["506","547","559"]:
                attempts = 0
                os.system(f'sleep=$(grep "Sleep" {params.batch_script} | cut -d " " -f 2) && sed "s/setSleepInterval $sleep/setSleepInterval 5000/g" {params.batch_script} > {params.batch_temp} && echo "exit" >> {params.batch_temp}')
                while height in ["506","547","559"] and attempts < 3:
                    os.system(f'echo "Snapshot may be truncated. Current snapshot height is {height}. Rerunning IGV batch script {params.batch_script} with increased sleep interval.\n" >> {log.stdout}')
                    attempts += 1
                    os.system(f'echo "IGV ATTEMPT #{attempts}:" >> {log.stdout}')
                    os.system(f'echo "IGV ATTEMPT #{attempts}:" >> {log.stderr}')
                    os.system(f'maxtime=$(($(wc -l < {params.batch_temp}) * 60 + 15)) ; timeout --foreground $maxtime xvfb-run -s "-screen 0 1980x1020x24" {params.server_number} {params.server_args} {params.igv} -b {params.batch_temp} >> {log.stdout} 2>> {log.stderr}')
                    height = str(subprocess.check_output(f"identify -format '%h' {input.snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                if height in ["547","559"]:
                    kurtosis, skewness = [float(value.split(": ")[1]) for value in str(subprocess.check_output(f'identify -verbose {input.snapshot} | grep -E "kurtosis|skewness" | tail -n 2', shell=True)).replace("\\n'","").split("\\n      ")]
                    blank_kurtosis = {"547": 18.5, "559": 18.2}
                    blank_skew = -4
                    if kurtosis > blank_kurtosis[height] and skewness < blank_skew:
                        os.system(f'echo "Snapshot may be blank. Current values:\nHeight:{height}, kurtosis: {str(kurtosis)}, skewness: {str(skewness)}\nSnapshots with height of 547, kurtosis greater than 18.5, and skewness less than 4 are likely blank, snapshots with height of 559, kurtosis greater than 18.2, and skewness less than 4 may be blank. Blank snapshots may be due to errors reading BAM file headers, Java address bind errors or other errors occurring during IGV run." >> {log.stdout}')
                        success = False
                if height == "506":
                    os.system(f'echo "Snapshot height is {height} and may still be truncated or improperly loaded. Check snapshot {input.snapshot}" >> {log.stdout}')
                    success = False
            if width == "640":
                attempts = 0
                os.system(f'echo "Snapshot appears to be in incorrect dimensions. Current width is {width} and might be due to xvfb-run unable to connect to current server {params.server_number} due to a server lock. Attempting to run on different server number..." >> {log.stdout}')
                if params.server_number == "--auto-servernum" or int(params.server_number.replace("-n ","")) >= 99:
                    new_server = 1
                else:
                    new_server = int(params.server_number.replace("-n ","")) + 1
                while width == "640" and attempts < 5:
                    os.system(f'echo "Attempting with server number {new_server}..." >> {log.stdout}')
                    os.system(f'echo "Attempting with server number {new_server}..." >> {log.stderr}')
                    os.system(f'maxtime=$(($(wc -l < {params.merged_batch}) * 60 + 15)) ; timeout --foreground $maxtime xvfb-run -s "-screen 0 1980x1020x24" -n {str(new_server)} {params.server_args} {params.igv} -b {params.merged_batch} >> {log.stdout} 2>> {log.stderr}')
                    width = str(subprocess.check_output(f"identify -format '%w' {input.snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                    new_server += 1
                if width == "640":
                    os.system(f'echo "Snapshot still appears to be in improper dimensions. Double check xvfb-run parameters." >> {log.stdout}')
                    success = False
            if success == True:
                os.system(f'touch {output.snapshot_qc}')
            if success == False:
                outline = "\t".join([wildcards.tumour_id, wildcards.seq_type, wildcards.genome_build, wildcards.gene, wildcards.chromosome, wildcards.start_position, wildcards.preset, input.snapshot])
                with open(params.failed_summary, "a") as handle:
                    handle.write(outline + "\n")

    # Symlinks the final output files into the module results directory (under '99-outputs/')
    rule _igv_symlink_snapshot:
        input:
            snapshot = CFG["dirs"]["snapshots"] + "{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".png",
            snapshot_qc = str(rules._igv_quality_control.output.snapshot_qc)
        output:
            snapshot = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}/{chromosome}:{start_position}--{gene}--{tumour_id}" + SUFFIX + ".png"
        threads:
            CFG["threads"]["_igv_symlink_snapshot"]
        run:
            op.relative_symlink(input.snapshot, output.snapshot)

    # Return a list of all snapshots that were taken during IGV for each specific tumour_id, tumour_seq_type, and tumour_genome_build combination
    def _symlink_snapshot(wildcards):
        CFG = config["lcr-modules"]["igv"]
        checkpoint_outputs = checkpoints._igv_run.get(**wildcards).output.complete

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.tumour_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)

        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf), 
            zip, 
            seq_type=wildcards.seq_type, 
            genome_build=wildcards.genome_build, 
            tumour_id=wildcards.tumour_id, 
            normal_sample_id=normal_sample_id, 
            pair_status=pair_status
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
                            tumour_id = maf_table["Tumor_Sample_Barcode"],
                            allow_missing = True
                        ),
                        pair_status_directory = PAIR_STATUS_DICT[pair_status.item()],
                        preset = wildcards.preset
                    )
        else:
            return []

    def _quality_control(wildcards):
        CFG = config["lcr-modules"]["igv"]

        checkpoint_outputs = checkpoints._igv_run.get(**wildcards).output.complete

        this_sample = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.tumour_id, tumour_seq_type = wildcards.seq_type, tumour_genome_build = wildcards.genome_build)

        normal_sample_id = this_sample["normal_sample_id"]
        pair_status = this_sample["pair_status"]

        maf = expand(
            str(rules._igv_filter_maf.output.maf),
            zip,
            seq_type = wildcards.seq_type,
            genome_build = wildcards.genome_build,
            tumour_id = wildcards.tumour_id,
            normal_sample_id = normal_sample_id,
            pair_status = pair_status
        )

        if os.path.exists(maf[0]):
            maf_table = pd.read_table(maf[0], comment="#", sep="\t")

            return expand(
                        expand(
                            str(rules._igv_quality_control.output.snapshot_qc),
                            zip,
                            seq_type = maf_table["seq_type"],
                            genome_build = maf_table["genome_build"],
                            chromosome = maf_table["chr_std"],
                            start_position = maf_table["Start_Position"],
                            gene = maf_table["Hugo_Symbol"],
                            tumour_id = maf_table["Tumor_Sample_Barcode"],
                            allow_missing = True
                        ),
                        pair_status_directory = PAIR_STATUS_DICT[pair_status.item()],
                        preset = wildcards.preset
                    )

    # Check that snapshots have been symlinked and quality controlled
    rule _igv_check_snapshots:
        input:
            igv_completed = str(rules._igv_run.output.complete),
            snapshots = _symlink_snapshot,
            quality_control = _quality_control
        output:
            snapshots = CFG["dirs"]["outputs"] + "completed/{preset}/{seq_type}--{genome_build}--{tumour_id}.completed"
        shell:
            "touch {output.snapshots}"

if CFG["estimate_only"] is True and CFG["identify_failed_snaps"] is False: 

    rule _igv_touch_summary:
        output:
            snapshot_summary = CFG["dirs"]["outputs"] + "snapshot_estimates/snapshot_summary.txt",
            snapshot_estimate = CFG["dirs"]["outputs"] + "snapshot_estimates/snapshot_estimate.txt",
            summary_ready = temp(CFG["dirs"]["outputs"] + "snapshot_estimates/touch_summary.completed")
        run:
            header = "\t".join(["sample_id","seq_type","genome_build","gene","chromosome","position", "igv_preset"])
            with open(output.snapshot_summary,"w") as handle:
                handle.write(header + "\n")
            with open(output.snapshot_estimate, "w") as handle:
                handle.write(header + "\n")
            ready = open(output.summary_ready, "w")
            ready.close()

    rule _igv_estimate_snapshots:
        input:
            maf = str(rules._igv_filter_maf.output.maf),
            summary_file_ready = str(rules._igv_touch_summary.output.summary_ready)
        output:
            estimate_finished = temp(CFG["dirs"]["batch_scripts"] + "estimate_batch_scripts/{seq_type}--{genome_build}/{preset}/{tumour_id}--{normal_sample_id}--{pair_status}.temp")
        params:
            snapshot_summary = str(rules._igv_touch_summary.output.snapshot_summary),
            snapshot_estimate = str(rules._igv_touch_summary.output.snapshot_estimate)
        threads: (workflow.cores)
        run:
            CFG = config["lcr-modules"]["igv"]
            maf_table = pd.read_table(input.maf, sep="\t", comment="#")

            seq_type = wildcards.seq_type
            genome_build = wildcards.genome_build
            tumour_id = wildcards.tumour_id
            preset = wildcards.preset

            snapshot_summary = open(params.snapshot_summary, "a")
            snapshot_estimate = open(params.snapshot_estimate, "a")

            for index, row in maf_table.iterrows():
                gene = row["Hugo_Symbol"]
                chromosome = row["chr_std"]
                position = str(row["Start_Position"])

                dispatch_path = CFG["dirs"]["batch_scripts"] + f"dispatched_batch_scripts/{seq_type}--{genome_build}/{preset}/{chromosome}:{position}--{gene}--{tumour_id}" + SUFFIX + ".batch"

                outline = "\t".join([tumour_id, seq_type, genome_build, gene, chromosome, position, preset])
                snapshot_summary.write(outline + "\n")

                if not os.path.exists(dispatch_path):
                    snapshot_estimate.write(outline + "\n")
            
            finished = open(output.estimate_finished, "w")
            finished.close()

if CFG["identify_failed_snaps"] is True and CFG["estimate_only"] is False:
    rule _igv_touch_failed:
        input:
            filter_maf = expand(str(rules._igv_filter_maf.output.maf), zip, seq_type = CFG["runs"]["tumour_seq_type"], tumour_id = CFG["runs"]["tumour_sample_id"], genome_build=CFG["runs"]["tumour_genome_build"], normal_sample_id=CFG["runs"]["normal_sample_id"], pair_status=CFG["runs"]["pair_status"])
        output:
            failed_summary = CFG["dirs"]["outputs"] + "snapshot_estimates/failed_summary.txt",
            failed_ready = temp(CFG["dirs"]["outputs"] + "snapshot_estimates/touch_failed.completed")
        run:
            header = "\t".join(["sample_id","seq_type","genome_build","gene","chromosome","position","preset","snapshot_path"])
            with open(output.failed_summary, "w") as handle:
                handle.write(header + "\n")
            ready = open(output.failed_ready, "w")
            ready.close()

    rule _igv_find_failed:
        input:
            maf = str(rules._igv_filter_maf.output.maf),
            failed_ready = str(rules._igv_touch_failed.output.failed_ready)
        output:
            failed_finished = temp(CFG["dirs"]["batch_scripts"] + "estimate_failed_scripts/{seq_type}--{genome_build}/{preset}/{tumour_id}--{normal_sample_id}--{pair_status}.temp")
        params:
            failed_summary = str(rules._igv_touch_failed.output.failed_summary)
        threads: (workflow.cores)
        run:
            import subprocess
            CFG = config["lcr-modules"]["igv"]
            maf_table = pd.read_table(input.maf, sep="\t", comment="#")

            seq_type = wildcards.seq_type
            genome_build = wildcards.genome_build
            tumour_id = wildcards.tumour_id
            preset = wildcards.preset
            pair_status_directory = PAIR_STATUS_DICT[wildcards.pair_status]

            for index, row in maf_table.iterrows():
                gene = row["Hugo_Symbol"]
                chromosome = row["chr_std"]
                position = str(row["Start_Position"])

                snapshot = CFG["dirs"]["snapshots"] + f"{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}/{chromosome}:{position}--{gene}--{tumour_id}" + SUFFIX + ".png"
                snapshot_symlink = CFG["dirs"]["outputs"] + f"{seq_type}--{genome_build}/{pair_status_directory}/{preset}/{chromosome}/{chromosome}:{position}--{gene}--{tumour_id}" + SUFFIX + ".png"

                success = True
                if not os.path.exists(snapshot):
                    print(f"{snapshot} doesn't exist yet, skipping...")
                if os.path.exists(snapshot):
                    if not os.path.exists(snapshot_symlink):
                        success = False
                    if success is True:
                        try:
                            height = str(subprocess.check_output(f"identify -format '%h' {snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                            width = str(subprocess.check_output(f"identify -format '%w' {snapshot}", shell=True)).split("'")[1].split("\\n")[0]
                            print(f"Height is {height} for {tumour_id} {gene} {chromosome}:{position}")
                        except:
                            success = False
                        if success is True and height in ["547","559"]:
                            kurtosis, skewness = [float(value.split(": ")[1]) for value in str(subprocess.check_output(f'identify -verbose {snapshot} | grep -E "kurtosis|skewness" | tail -n 2', shell=True)).replace("\\n'","").split("\\n      ")]
                            print(f"Kurtosis is {kurtosis} for {tumour_id} {gene} {chromosome}:{position}")
                            blank_kurtosis = {"547": 18.5, "559": 18.2}
                            blank_skew = -4
                            if kurtosis > blank_kurtosis[height] and skewness < blank_skew:
                                success = False
                print(f"Success value is {success} for {tumour_id} {gene} {chromosome}:{position}")
                if success is False:
                    with open(params.failed_summary, "a") as handle:
                        print("Writing line to failed file")
                        outline = "\t".join([tumour_id, seq_type, genome_build, gene, chromosome, position, preset, snapshot])
                        handle.write(outline + "\n")
                    
            finished = open(output.failed_finished, "w")
            finished.close()

# Generates the target sentinels for each run, which generate the symlinks
if CFG["estimate_only"] is False and CFG["identify_failed_snaps"] is False:
    rule _igv_all:
        input:
            expand(
                expand(
                    [
                        str(rules._igv_run.output.complete), 
                        str(rules._igv_check_snapshots.output.snapshots)
                    ], 
                    zip, 
                    tumour_id=CFG["runs"]["tumour_sample_id"], 
                    seq_type=CFG["runs"]["tumour_seq_type"], 
                    genome_build=CFG["runs"]["tumour_genome_build"],
                    allow_missing=True
                ),
                preset=CFG["options"]["igv_presets"]
            )

if CFG["estimate_only"] is True and CFG["identify_failed_snaps"] is False:
    rule _igv_all:
        input:
            expand(
                expand(
                    str(rules._igv_estimate_snapshots.output.estimate_finished),
                    zip,
                    seq_type=CFG["runs"]["tumour_seq_type"],
                    genome_build=CFG["runs"]["tumour_genome_build"],
                    tumour_id=CFG["runs"]["tumour_sample_id"],
                    normal_sample_id=CFG["runs"]["normal_sample_id"],
                    pair_status=CFG["runs"]["pair_status"],
                    allow_missing=True
                ),
                preset=CFG["options"]["igv_presets"]
            )

if CFG["identify_failed_snaps"] is True and CFG["estimate_only"] is False:
    rule _igv_all:
        input:
            expand(
                expand(
                    str(rules._igv_find_failed.output.failed_finished),
                    zip,
                    seq_type=CFG["runs"]["tumour_seq_type"],
                    genome_build=CFG["runs"]["tumour_genome_build"],
                    tumour_id=CFG["runs"]["tumour_sample_id"],
                    normal_sample_id=CFG["runs"]["normal_sample_id"],
                    pair_status=CFG["runs"]["pair_status"],
                    allow_missing=True
                ),
                preset=CFG["options"]["igv_presets"]
            )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
