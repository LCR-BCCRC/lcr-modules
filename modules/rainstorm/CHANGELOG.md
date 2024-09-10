# Changelog

All notable changes to the `rainstorm` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-05-06

This release was authored by Kostiantyn Dreval.

- At the first implementation, the module runs both rainstorm and doppler analyses on a cohort-wide maf file.
- Installation of github repo is achieved through download of tar file with specific version tag, rather than
  repo cloning to ensure tool version consistency.

## [1.1] - 2024-07-30

This update was authored by Sierra Gillis.

- Added a preprocessing step that uses subsetting categories to subset the input maf to samples of interest. This is implemented as a checkpoint that creates files with md5sum in their name. The md5sum is of the sample IDs in the `sample_set` to allow tracking as this set of samples changes.
- Launch date is also used in file names to make identifying runs for different samples *within* the same `sample_set` easier at a glance. This way users do not need to only rely on the md5sum to identify their run. The sample IDs are also output to a file named with the md5sum so users can double check the correspondence between the two.
- There is an aggregate rule to trigger the checkpoint and remove the wildcard before `rule all`.
- Preprocessing step was updated to be agnostic of the input metadata samples dataframe. It will check the input file name and the subsetting categories to ensure that capture data is not used, as this module only runs on genome data.
- The workflow was modified to run on genome builds with and without "chr" prefixes, and names files accordingly.

## [1.2] - 2024-08-22

This update was authored by Sierra Gillis.

- Added an output directory for the prepare step
- Added more wildcards to the input maf symlink path for better knowing which was used for which run
- Was tested with an update made to `lcr-scripts/generate_smg_inputs` that uses threads to read in maf files faster, and only reads in relevant columns
