# Changelog

All notable changes to the `salmon` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-07-16

This release was authored by Helena Winata.

- Need link to download rna fasta file, `GRCh38_latest_rna.fna`, to generate salmon index

## [1.1] - 2020-10-30
This release was authored by Laura Hilton

- Optional job grouping is included in the snakefile `salmon_grouped.smk`. This enables rapid temp fastq deletion when used in conjunction with the `bam2fastq/1.2/bam2fastq_grouped.smk` module. 
- Utilizes resource unpacking for flexibility in specifying resource requirements. 