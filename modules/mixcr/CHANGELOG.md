# Changelog

All notable changes to the `mixcr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2] - 2022-06-6

This release was updated by Manuela Cruz.

- Removes dependency on `{genome_build}`wildcard.
- Adds requirement to specify desired receptor results in config.
- Explicitly defines output files for desired receptors as `mixcr.{sample_id}.clonotypes.{chain}.txt`. This allows `mixcr` to rerun if necessary to create desired receptors files.
- Adds option to run `igblastn` on `mixcr` results.
- Checks that `mixcr` run parameters are compatible with IgBLASTN requirements.
- Modifies run parameters if necessary to obtain results for desired receptors.
- Includes scripts for parsing `mixcr` and `igblastn` results and merging results into file `mixcr.{sample_id}.clonotypes.{chain}.status.txt`

## [1.1] - 2020-12-21

This release was updated by Laura Hilton. 

- The user now specifies where `mixcr` should be installed to prevent unwanted installation into the `lcr-modules` repository. 
- Utilizes resource unpacking. 
- Added a conda env to ensure Java > 8 is used. 

## [1.0] - 2020-06-11

This release was authored by Anita dos Santos.

- The `mixcr` module can process RNA-seq or non-targeted genomic raw data
- Default options are used for all optional arguments.
- Final outputs are named `mixcr.{sample_id}.clonotypes.ALL.txt` and `mixcr.{sample_id}.report` and are symlinked to `99-outputs`
- A `mixcr.{sample_id}.clonotypes.ALL.txt` file will be created regardless of the success of `mixcr` to circumvent cases where no clonotypes are found and no such file is created.
- The latest release of MiXCR is installed from the milaboratory MiXCR github as a jar. A conda environment for `mixcr` was forecast to be released in July 2020 but is not yet available as of Dec. 2020.