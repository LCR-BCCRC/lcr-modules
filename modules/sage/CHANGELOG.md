# Changelog

All notable changes to the `sage` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2022-04-26

This release was authored by Kostiantyn Dreval.

- The version 1.1 of the module uses version 3.0.1 of the SAGE variant caller. It was incompatible with the version 2.6 used before because argument names (`assembly` -> `ref_genome_version`; `chr` -> `specific_chr`) and the expected inputs of previous arguments (hg19/hg38  for ref_genome_version became 37/38; chromosome names for specific_chr now need to be specified through `;` instead of `,` ) got changed.
- The new required argument `ensembl_data_dir` was added, and in the module it is implemented as part of the rule `_download_sage_references`
- The advantage of this version is that it can now handle split reads where the length of deletion is longer than the read itself. Previous versio of the SAGE was erroring out on such instances, and the current version converts this to warning and skips such reads. For details see this issue: https://github.com/hartwigmedical/hmftools/issues/220

## [1.0] - 2020-12-29

This release was authored by Kostiantyn Dreval.

- Reference files are symlinked to `00-inputs` from references directory to ensure all required indexes and dictionaries are created if starting on fresh directory.
- SAGE fails to write to .gz files directly, therefore outputs are marked as `temp()` to ensure they are deleted after bgzip.
- `main_chromosomes.txt` from reference_files workflow is used to restrict SAGE only to standard chromosomes. Chromosomes listed in this file are returned as list to be included with SAGE `-chr` flag.
- module outputs separate VCF files for indels, snvs, and combined (indels+snsvs+mnvs). This can be modified in later releases.
- additional flags for run can be specified in `CFG["options"]["sage_run"]` without need to modify module snakefile.
