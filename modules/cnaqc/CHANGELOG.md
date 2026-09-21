# Changelog

All notable changes to the `cnaqc` module will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-05-27

### Added
- Initial implementation of the CNAqc module
- Takes Battenberg subclones (filled) and cellularity_ploidy as CNA inputs
- Takes somatic SNV MAF files (e.g. from slms-3/vcf2maf) as mutation inputs
- Runs CNAqc peak analysis and computes CNAqc QC score
- Outputs: per-sample QC PDF plots and metrics TSV
