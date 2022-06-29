# Changelog

All notable changes to the `oncodriveclustl` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-10-15

This release was authored by Manuela Cruz.

- The `oncodriveclustl` module takes in a mutations/MAF file and annotations/regions file to identify significantly clustered signals in genomic regions.
- MAF files are combined and converted into the Oncodrive format with an R script (samples with multiple seq_types and duplicate mutation calls are also handled by the R script).
- Running OncodriveCLUSTL:
    - Regions files must be specified in the config for each desired genome build.
    - A `local path` must be specified in the config. This determines where the `genome reference` and `genomic scores` directories are downloaded from `bbglab`. If this path is not specified, the default download path is `~/.bgdata/`. These are large files (~ 26G) so it is recommended to have a designated location for them.
    - All default parameters will be used unless explicitly specified in the config `CFG["options"]["clustl"]`
- Three output files are produced in the `01-clustl/{projection}/{sample_set}` directory: `clusters_results.tsv`, `elements_results.txt` and `quantile_quantile_plot.png`. They are symlinked to 99-outputs as `{sample_set}_clusters_results.tsv`, `{sample_set}_elements_results.txt`, and `{sample_set}_quantile_quantile_plot.png`.
