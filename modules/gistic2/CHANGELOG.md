# Changelog

All notable changes to the `gistic2` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-06-06

This release was authored by Sierra Gillis.

- All preprocess of the seg file is performed by one step.
Includes keeping genomic data over capture data where available, filtering out overlapping
segments, filtering out data for non-canonical chromosomes

## [1.1] - 2023-10-26

This update was authored by Sierra Gillis.

- Preprocessing step was converted to a checkpoint that creates files with md5sum in their name, used as a wildcard. The md5sum is of the sample IDs in the `case_set`
- There is an aggregate rule to trigger the checkpoint and remove the wildcard before `rule all`
- R Script in the preprocessing step no longer uses `argparase`, instead it uses snakemake objects. Particularly, the metadata for getting the sample_IDs in the `case_set` is the dataframe stored in `CFG["samples"]` subset to the appropriate columns.\
- Launch date is also used in file names to make identifying runs for different samples *within* the same `case_set` easier at a glance. This way users do not need to only rely on the md5sum to identify their run. The sample IDs are also output to a file named with the md5sum so users can double check the correspondence between the two.