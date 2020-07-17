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
