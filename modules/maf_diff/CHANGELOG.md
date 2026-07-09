# Changelog

All notable changes to the `maf_diff` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-09

### Added

- Initial implementation of `maf_diff` module
- Compares MAF files from two variant callers per tumour/normal pair
- Outputs caller-specific MAFs for variants unique to each caller
- Outputs a summary TSV with counts of shared, caller1-only, and caller2-only variants
- Matching is performed on configurable key columns (default: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
