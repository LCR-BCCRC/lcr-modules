# Changelog

All notable changes to the `dnds` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-08-02

This release was authored by Kostia Dreval.

- The version 1 uses the reference files supplied directly with dNdS.
- Currently, only the grch37-based maf files are supported.
- The snakefiles relies on the lcr-script `generate_smg_inputs` to produce the input maf.
- The implementation is similar to the other SMG lcr-modules, e.g. mutsig.

## [1.1] - 2024-01-

This update was authored by Sierra Gillis.

- Preprocessing maf step was updated to be agnostic of the input metadata samples dataframe. The values to subset this metadata by in order to get the IDs of samples in the sample set will be done through the level3_subsetting_categories.tsv table.
- The sample set outputs are tracked by launch-date of the snakemake and md5sum of the sample ids as well as sample set name. This implementation follows that in gistic2/1.1/

## [1.2] - 2024-08-22

This update was authored by Sierra Gillis.

- Added an output directory for the prepare step
- Added more wildcards to the input maf symlink path for better knowing which was used for which run
- Was tested with an update made to `lcr-scripts/generate_smg_inputs` that uses threads to read in maf files faster