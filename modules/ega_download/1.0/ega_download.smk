#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
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
    subdirectories = ["inputs", "ega_download", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _ega_input_csv,
    _ega_output_files,
    _ega_all


##### RULES #####

# Confirm the output file type is correctly specified.
out_file_type = CFG["out_file_type"]
possible_types = ["bam", "cram", "fastq", "fq.gz", "fastq.gz"]
possible_type_string = ", ".join(possible_types)
assert out_file_type in possible_types, (
    f"Unrecognized out_file_type {out_file_type}. \n"
    f"Please specify one of {possible_type_string}. "
)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _ega_input_csv:
    input:
        csv = ancient(CFG["inputs"]["sample_csv"])
    output:
        csv = CFG["dirs"]["inputs"] + "{study_id}.csv"
    run:
        op.absolute_symlink(input.csv, output.csv)

# Download the file from EGA
rule _ega_get_ega_file:
    input:
        sample_table = str(rules._ega_input_csv.output.csv)
    output:
        ega_file = CFG["dirs"]["ega_download"] + "{seq_type}/{file_format}/{study_id}/{egaf}/{egas}--{egad}--{egan}--{egaf}--{file_name}.{file_format}"
    log:
        stdout = CFG["logs"]["ega_download"] + "{seq_type}/{file_format}/{study_id}/{egaf}/{egas}--{egad}--{egan}--{egaf}--{file_name}.{file_format}_download.stdout.log",
        stderr = CFG["logs"]["ega_download"] + "{seq_type}/{file_format}/{study_id}/{egaf}/{egas}--{egad}--{egan}--{egaf}--{file_name}.{file_format}_download.stderr.log"
    params:
        credentials = CFG["credentials_file"],
        additional_args = CFG["options"]["additional_args"]
    conda:
        CFG["conda_envs"]["pyega3"]
    threads:
        CFG["threads"]["ega_file_download"]
    resources:
        **CFG["resources"]["ega_file_download"]
    wildcard_constraints:
        file_format = "|".join(CFG["samples"]["file_format"].tolist())
    shell:
        op.as_one_line("""
        pyega3
        {params.additional_args}
        -cf {params.credentials}
        fetch {wildcards.egaf}
        --output-dir $(dirname "$(dirname {output.ega_file})")
        --max-retries -1
        --retry-wait 30
        > {log.stdout}
        2> {log.stderr}
        &&
        mv $(dirname {output.ega_file})/{wildcards.file_name}.{wildcards.file_format} {output.ega_file}
        """)

# function to get ega_file files that drops EGA IDs from the name of the final
# file and uses sample_id instead
def get_ega_file (wildcards):
    CFG = config["lcr-modules"]["ega_download"]
    tbl = CFG["samples"]
    this_sample = tbl[(tbl.seq_type == wildcards.seq_type) & (tbl.sample_id == wildcards.sample_id) & (tbl.file_format == wildcards.file_format)]
    if (this_sample.shape[0]>0):
        this_sample = this_sample.iloc[int(wildcards.read) - 1]
    this_egas, this_egad, this_egan, this_egaf, this_name = this_sample['EGAS'], this_sample['EGAD'], this_sample['EGAN'], this_sample['EGAF'], this_sample['file_name']
    this_file = expand(
            str(
                CFG["dirs"]["ega_download"] + "{{seq_type}}/{{file_format}}/{{study_id}}/{egaf}/{egas}--{egad}--{egan}--{egaf}--{{file_name}}.{{file_format}}"
            ),
            egas = this_egas,
            egad = this_egad,
            egan = this_egan,
            egaf = this_egaf
        )[0]
    this_file = str(this_file).replace("{file_name}", this_name)
    return(this_file)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ega_output_files:
    input:
        ega_file = get_ega_file
    output:
        ega_file = CFG["dirs"]["outputs"] + "{seq_type}/{study_id}/{sample_id}{read}.{file_format}"
    wildcard_constraints:
        file_format = "|".join(CFG["samples"]["file_format"].tolist())
    run:
        op.relative_symlink(input.ega_file, output.ega_file, in_module=True)

# Generates the target files for each run, which generate the symlinks
rule _ega_all:
    input:
        expand(
            [
                str(rules._ega_output_files.output.ega_file)
            ],
            zip,  # Run expand() with zip(), not product()
            sample_id=CFG["samples"]["sample_id"],
            file_format=CFG["samples"]["file_format"],
            seq_type=CFG["samples"]["seq_type"],
            study_id=[CFG["study_id"]]*len(CFG["samples"]["file_name"]),
            read = [""]
        ) if out_file_type in ["cram", "bam"] else expand(
            expand(
                [
                    str(rules._ega_output_files.output.ega_file)
                ],
                zip,  # Run expand() with zip(), not product()
                sample_id=CFG["samples"]["sample_id"],
                file_format=CFG["samples"]["file_format"],
                seq_type=CFG["samples"]["seq_type"],
                study_id=[CFG["study_id"]]*len(CFG["samples"]["file_name"]),
                allow_missing = True
            ),
        read = ["1", "2"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
