# Changelog

All notable changes to the `oncodrivefml` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-10-25

This release was authored by Manuela Cruz.

- Input MAFs for OncodriveFML are created by the lcr-scripts `generate_smg_inputs`. Desired cohorts, pathologies, genome builds, etc. can be defined within a subsetting file and provided in the config at ["inputs"]["subsetting_categories"]

- Additional blacklists can be provided to filter mutations in ["maf_processing"]["blacklists"]

- OncodriveFML depends on regions files. grch37 regions files are available through the `reference_files` workflow. Available regions are cds, lincrna, 5utr, 3utr and promoter_splice_sites. The hg38 cds regions file is available within the module at `etc/hg38.cds.tsv.gz`. The path to the reference_files directory must be provided within the config at ["reference_files_directory"] so that the global variable `BGDATA_LOCAL` can be set.

- OncodriveFML also depends on the hg19 or hg38 genome reference files, which are available through the `reference_files` workflow as well.

- The `BGDATA_LOCAL` variable is required by the `bgdata` package within OncodriveFML. If this variable is not set, the default path is `~/.bgdata/` and will result in the genome reference and regions files being downloaded at this location. It is important to set this to the correct reference files path to prevent unnecessary downloads.

- OncodriveFML depends on CADD scores to score the deleteriousness of each mutation. The grch37 CADD 1.0 fml-formatted scores are provided by the `bgdata` package. This download will occur in the `01-cadd` scratch directory as the file is 19G. If you wish to use scores for a different build, you can provide a path to the score file in the config at ["options"]["score_path"]. The directory path must also be provided in a separate config file, detailed in the next bullet point. More information on CADD and downloading score packs here: https://cadd.gs.washington.edu/

- OncodriveFML depends on an additional config file to change analysis behaviour. The default configuration file is provided in the module directory at `etc/oncodrivefml.conf`. You must provide the path to the directory of the scores file at the `file =` option. More information on OncodriveFML config parameters here: https://oncodrivefml.readthedocs.io/en/latest/configuration.html
