#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import pandas as pd
import json
import snakemake.remote.EGA as EGA
import os
import stat

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
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["ega_download"]`
CFG = op.setup_module(
    name = "ega_download",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "ega_download", "outputs"],
)

# Check if EGA credentials file is only accessible by user or user group. If someone outside user group can read/write this file, stop module from the execution
permissions = (oct(stat.S_IMODE(os.lstat(CFG["credentials_file"]).st_mode)))[-1:]
if (int(permissions) != 0):
    sys.exit("The EGA credentials file is readable/writebale not only by the owner or user group. Please ensure that EGA credentials file can only be accessed by the owner or its user group.")

# Opening JSON file with credentials
EGA_CREDENTIALS_FILE = open(CFG["credentials_file"],)
# Returns credentials as a dictionary
EGA_CREDENTIALS = json.load(EGA_CREDENTIALS_FILE )
# set env variables to connect to EGA
os.environ["EGA_USERNAME"] = EGA_CREDENTIALS['username']
os.environ["EGA_PASSWORD"] = EGA_CREDENTIALS['password']
os.environ["EGA_CLIENT_SECRET"] = "AMenuDLjVdVo4BSwi0QD54LL6NeVDEZRzEQUJ7hJOM3g4imDZBHHX0hNfKHPeQIGkskhtCmqAJtt_jm7EKq-rWw"
os.environ["EGA_CLIENT_ID"] = "f20cd2d3-682a-4568-a53e-4262ef54c8f4"
# Closing file
EGA_CREDENTIALS_FILE.close()


# Defining global variables and doing global setup to connect to EGA
ega = EGA.RemoteProvider()
# The file listing samples
EGA_TARGET_SAMPLES = pd.read_csv (str(CFG["inputs"]["sample_csv"]).replace("{study_id}", CFG["study_id"]))
# The file with metadata
EGA_TARGET_SAMPLES_METADATA = pd.read_csv (str(CFG["inputs"]["sample_metadata"]).replace("{study_id}", CFG["study_id"]))
# Merge them together, making sure the sample table is sorted
EGA_MASTER_SAMPLE_TABLE = EGA_TARGET_SAMPLES.merge(EGA_TARGET_SAMPLES_METADATA, on="file_id", how="right")
EGA_MASTER_SAMPLE_TABLE = EGA_MASTER_SAMPLE_TABLE.sort_values(by=['file_id'], ascending=True)
tbl = EGA_MASTER_SAMPLE_TABLE

# Define rules to be run locally when using a compute cluster
localrules:
    _ega_download_input_csv,
    _ega_download_get_bam,
    _ega_download_index_bam,
    _ega_download_output_files,
    _ega_download_export_sample_table,
    _ega_download_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ega_download_input_csv:
    input:
        csv = ancient(CFG["inputs"]["sample_csv"]),
        metadata = ancient(CFG["inputs"]["sample_metadata"])
    output:
        csv = CFG["dirs"]["inputs"] + "csv/{seq_type}--{genome_build}/{study_id}.csv",
        metadata = CFG["dirs"]["inputs"] + "metadata/{seq_type}--{genome_build}/{study_id}.metadata.csv",
    run:
        op.absolute_symlink(input.csv, output.csv)
        op.absolute_symlink(input.metadata, output.metadata)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _ega_download_get_bam:
    input:
        sample_table = str(rules._ega_download_input_csv.output.csv),
        #bam = ega.remote(str("ega/{study_id}/{sample_id}.bam"))
        bam = str("{study_id}/{file_id}.bam")
    output:
        bam = CFG["dirs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}.bam"
    log:
        stdout = CFG["logs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}_download.stdout.log",
        stderr = CFG["logs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}_download.stderr.log"
    threads:
        CFG["threads"]["bam_index"]
    resources:
        **CFG["resources"]["bam_download"]
    shell:
        "cp {input.bam} {output.bam} > {log.stdout} 2> {log.stderr}"


# Example variant filtering rule (single-threaded; can be run on cluster head node)
rule _ega_download_index_bam:
    input:
        bam = str(rules._ega_download_get_bam.output.bam)
    output:
        bai = CFG["dirs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}.bam.bai"
    log:
        stdout = CFG["logs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}_index.stdout.log",
        stderr = CFG["logs"]["ega_download"] + "bam/{seq_type}--{genome_build}/{study_id}/{file_id}/{sample_id}_index.stderr.log"
    params:
        opts = CFG["options"]["bam_index"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["bam_index"]
    resources:
        **CFG["resources"]["bam_index"]
    shell:
        op.as_one_line("""
        samtools index
        {params.opts}
        -@ {threads}
        {input.bam}
        > {log.stdout}
        2> {log.stderr}
        """)


# function to get bam files that drops EGA ID from the name of the final bam
def get_file_bam (wildcards):
    CFG = config["lcr-modules"]["ega_download"]
    file_id = tbl[(tbl.sample_id == wildcards.sample_id) & (tbl.seq_type == wildcards.seq_type)]["file_id"]
    return (expand(str(CFG["dirs"]["ega_download"] + "bam/{{seq_type}}--{{genome_build}}/{{study_id}}/{file_id}/{{sample_id}}.bam"),
                    file_id = file_id))

# function to get bam files that drops EGA ID from the name of the index file
def get_file_bai (wildcards):
    CFG = config["lcr-modules"]["ega_download"]
    file_id = tbl[(tbl.sample_id == wildcards.sample_id) & (tbl.seq_type == wildcards.seq_type)]["file_id"]
    return (expand(str(CFG["dirs"]["ega_download"] + "bam/{{seq_type}}--{{genome_build}}/{{study_id}}/{file_id}/{{sample_id}}.bam.bai"),
                    file_id = file_id))

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ega_download_output_files:
    input:
        bam = get_file_bam,
        bai = get_file_bai
    output:
        bam = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{study_id}/{sample_id}.bam",
        bai = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{study_id}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam, in_module=True)
        op.relative_symlink(input.bai, output.bai, in_module=True)


# Write associated sample table as one of the module's outputs
rule _ega_download_export_sample_table:
    input:
        expand(
            [
                str(rules._ega_download_output_files.output.bam),
                str(rules._ega_download_output_files.output.bai)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=EGA_MASTER_SAMPLE_TABLE["seq_type"],
            genome_build=EGA_MASTER_SAMPLE_TABLE["genome_build"],
            study_id=[CFG["study_id"]]*len(EGA_TARGET_SAMPLES["file_id"]),
            sample_id=EGA_MASTER_SAMPLE_TABLE["sample_id"])
    output:
        sample_table = CFG["dirs"]["outputs"] + "metadata/{study_id}/metadata.tsv",
    run:
        samples_metadata = EGA_MASTER_SAMPLE_TABLE.rename(columns={'file_format': 'compression'})
        samples_metadata = samples_metadata.insert(0, 'bam_available', 'TRUE')
        samples_metadata.to_csv(output.sample_table, sep='\t', header=True, index=False)


# Generates the target files for each run, which generate the symlinks
rule _ega_download_all:
    input:
        expand(
            [
                str(rules._ega_download_output_files.output.bam),
                str(rules._ega_download_output_files.output.bai),
                str(rules._ega_download_export_sample_table.output.sample_table)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=EGA_MASTER_SAMPLE_TABLE["seq_type"],
            genome_build=EGA_MASTER_SAMPLE_TABLE["genome_build"],
            study_id=[CFG["study_id"]]*len(EGA_TARGET_SAMPLES["file_id"]),
            sample_id=EGA_MASTER_SAMPLE_TABLE["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
