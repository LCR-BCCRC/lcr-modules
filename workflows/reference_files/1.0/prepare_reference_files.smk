#!/usr/bin/env snakemake


##### REFERENCE_FILES SUBWORKFLOW #####


subworkflow reference_files:
    workdir:
        "/projects/bgrande/reference"
    snakefile:
        "reference_files.smk"
    configfile:
        "config/default.yaml"


##### REFERENCE_FILES TARGETS #####


# Must load configfile before including `reference_files.smk`
# Also used for wildcard values and tool versions
configfile: "config/default.yaml"


# Ask for all copies of all output files in `reference_files.smk`
rule all:
    input:
        reference_files(expand(
            [
                "genomes/{build}/genome_fasta/genome.fa",
                "genomes/{build}/genome_fasta/genome.fa.fai",
                "genomes/{build}/genome_fasta/main_chromosomes.txt",
                "genomes/{build}/genome_fasta/main_chromosomes.bed",
                "genomes/{build}/bwa_index/bwa-{bwa_version}/genome.fa",
                "genomes/{build}/annotations/gencode_annotation-{gencode_release}.gtf",
                "genomes/{build}/star_index/star-{star_version}/gencode-{gencode_release}/overhang-{star_overhang}",
            ],
            build=config["genome_builds"].keys(),
            bwa_version=config["tools"]["bwa"]["version"],
            gencode_release=config["wildcard_values"]["gencode_release"],
            star_version=config["tools"]["star"]["version"],
            star_overhang=config["wildcard_values"]["star_overhang"],
        ))
