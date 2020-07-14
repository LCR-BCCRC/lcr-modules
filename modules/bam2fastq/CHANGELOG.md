# Changelog

All notable changes to the `bam2fastq` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-29

This release was authored by Helena Winata.


- The `bam2fastq` convert bam files to paired end fastqs.
- When using in conjuction with a module that takes fastq inputs users have the option to remove fastq files after.
    - This is done via the `_bam2fastq_delete_fastq` and `_bam2fastq_delete_fastq_all`.
    - It follows the mechanism of a temp file where the file is deleted when a certain output is created.
    - Users need to specify the output in `config["lcr-modules"]["bam2fastq"]["outputs"]` to ensure that the fastq files are deleted only after a certain output is created.
