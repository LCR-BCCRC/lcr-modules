# Changelog

All notable changes to the `mutect2` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-08-15

This release was authored by Prasath Pararajalingam.
- Added rules to run tumours in no normal mode.
- Added parallelization by chromosome.
- gnomAD VCFs are chr-prefixed depending on the genome build. CVbio rules in `reference_files_header.smk` expect a single provider value for updating contig names between builds, however, the gnomAD providers change based on genome build. `get_cvbio_params()` and `hardlink_same_provider()` were updated to handle functions which can return a single value based on the genome version.
- Added reference rules to download the AF only gnomAD VCF file from [here](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/).
- Added reference rules to create .dict file within genome directory. The fasta index (.fai) and .dict file must be in the same directory.
- The BAM SM tag may be different from the BAM filename. Mutect2 requires the correct SM tag specification. `_mutect2_get_sm` was created to retrieve and output the SM tag into a file.
- Outputs compressed, PASS-filtered variants in final output directory. Unfiltered variants output into `filter` directory.
