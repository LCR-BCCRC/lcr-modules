# Changelog

All notable changes to the `mixcr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-12-08

This release was authored by Anita dos Santos with contributions from Laura Hilton.

- The `mixcr` module can process RNA-seq or non-targeted genomic raw data
- Default options are used for all optional arguments.
- Final outputs are named `mixcr.{sample_id}.clonotypes.ALL.txt` and `mixcr.{sample_id}.report` and are symlinked to `99-outputs`
- A `mixcr.{sample_id}.clonotypes.ALL.txt` file will be created regardless of the success of `mixcr` to circumvent cases where no clonotypes are found and no such file is created.
- The latest release of MiXCR is installed from the milaboratory MiXCR github as a jar. A conda environment for `mixcr` was forecast to be released in July 2020 but is not yet available as of Dec. 2020.