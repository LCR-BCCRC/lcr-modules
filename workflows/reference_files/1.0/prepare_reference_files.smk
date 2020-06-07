#!/usr/bin/env snakemake


##### SETUP #####


import os


##### REFERENCE_FILES MODULE #####


# Must load configfile before including `reference_files.smk`
# Also used for wildcard values and tool versions
configfile: "config/default.yaml"

# Check that the reference_directory is provided
assert "reference_directory" in config, (
    "Re-run snakemake with `--config reference_directory=/path/to/reference_files` in your command."
)

# Confirm that current working directory matched config value
original_dir = os.getcwd()
reference_dir = config["reference_directory"]
os.chdir(reference_dir)

# Include the `reference_files` module
include: os.path.join(original_dir, "reference_files.smk")


##### REFERENCE_FILES TARGETS #####


# Ask for all copies of all output files in `reference_files.smk`
rule all:
    input:
        expand(
            [
                rules.download_genome_gc.output.gc,
                rules.download_genome_dbsnp.output.dbsnp,
                rules.download_genome_gnomad.output.gnomad,
                rules.get_genome_fasta_download.output.fasta,
                rules.index_genome_fasta.output.fai,
                rules.get_main_chromosomes_download.output.txt,
                rules.get_main_chromosomes_download.output.bed,
                rules.create_bwa_index.output.prefix,
                rules.get_gencode_download.output.gtf,
                rules.create_star_index.output.index,
            ],
            genome_build=config["genome_builds"].keys(),
            bwa_version=config["tools"]["bwa"]["version"],
            gencode_release=config["wildcard_values"]["gencode_release"],
            star_version=config["tools"]["star"]["version"],
            star_overhang=config["wildcard_values"]["star_overhang"],
        )
