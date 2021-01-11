# Changelog

All notable changes to the `sage` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-12-29

This release was authored by Kostiantyn Dreval.

- Reference files are symlinked to `00-inputs` from references directory to ensure all required indexes and dictionaries are created if starting on fresh directory.
- SAGE fails to write to .gz files directly, therefore outputs are marked as `temp()` to ensure they are deleted after bgzip.
- `main_chromosomes.txt` from reference_files workflow is used to restrict SAGE only to standard chromosomes. Chromosomes listed in this file are returned as list to be included with SAGE `-chr` flag.
- module outputs separate VCF files for indels, snvs, and combined (indels+snsvs+mnvs). This can be modified in later releases.
- additional flags for run can be specified in `CFG["options"]["sage_run"]` without need to modify module snakefile.