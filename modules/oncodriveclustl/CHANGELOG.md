# Changelog

All notable changes to the `oncodriveclustl` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2024-02-07

This release was authored by Manuela Cruz.

- Input MAFs for OncodriveCLUSTL are created by the lcr-scripts `generate_smg_inputs`. Desired cohorts, pathologies, genome builds, etc. can be defined and provided in the config at ["inputs]["subsetting_categories"]

- Additional blacklists can be provided to filter mutations in ["maf_processing"]["blacklists"]

- OncodriveCLUSTL depends on regions files. grch37 regions files are available through the `reference_files` workflow. Available regions are cds, lincrna, 5utr, 3utr and 5utr. The hg38 cds region file is available within the module at `etc/hg38.cds.tsv.gz`. The path to the reference_files directory must be provided within the config at ["reference_files_directory"] so that the global variable `BGDATA_LOCAL` can be set.

- OncodriveCLUSTL also depends on hg19 or hg38 genome reference files, which are available through the `reference_files` workflow as well.

- The `BGDATA_LOCAL` variable is required by the `bgdata` package within OncodriveCLUSTL. If this variable is not set, the default path is `~/.bgdata/` and will result in the genome reference and regions files being downloaded at this location. It is important to set this to the correct reference files path to prevent unnecessary downloads.

- The script `etc/detailed_clusters.py` expands the returned clusters into genomic positions. Desired clusters are filtered by setting desired q-values of the elements using the config value ["q_values"] and desired p-values, minimum_samples, and minimum_score in the config at ["detailed_cluster_options"]