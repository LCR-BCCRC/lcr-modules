# Changelog

All notable changes to the `bam2fastq` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-29

This release was authored by Helena Winata.


- The `bam2fastq` convert bam files to paired end fastqs.
- When using in conjuction with a module that takes fastq inputs users have the option to make fastq files temporary.
    - This is done via by setting `CFG["temp_outputs"]` to `True`.
    - This instructs snakemake to use a version of the `_bam2fastq_run` rule where the outputs are temporary

## [1.1] - 2020-08-25

Tweaks by Ryan Morin

- We need to start handling CRAM-compressed bam files. This modification works on bam or cram files (regardless of their name).
- The new version now relies on the reference_files function because the genome fasta is needed for cram decompression
- tested on some very old RNA-seq data (DLBCL cell lines)

- The outputs are not separated by genome build so that it can be re-aligned using other genome references in later modules

On 2020-09-11, the new updates were introduced without version increase. Mainly, the genome build wildcard
was removed from the module using a different approach developed by the team. In addition, the fastq files
are now rule outputs, and therefore subsequent modules can recognize these files being produced.