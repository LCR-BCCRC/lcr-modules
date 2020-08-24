# Changelog

All notable changes to the `mutect2` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-08-15

This release was authored by Prasath Pararajalingam.

- Added reference rules to create .dict file within genome directory. The fasta index (.fai) and .dict file must be in the same directory.
- The BAM SM tag may be different from the BAM filename. Mutect2 requires the correct SM tag specification. `_mutect2_get_sm` was created to retrieve and output the SM tag into a file.
- Added reference rules to download the AF only gnomAD VCF file from [here](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2/).
- Outputs compressed, PASS-filtered variants in final output directory. Unfiltered variants output into `XX-filter` directory.
