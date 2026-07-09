# Changelog

All notable changes to the `igv_reports` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-09

### Added

- Initial implementation of `igv_reports` module wrapping igvteam/igv-reports v1.16.0
- Filters input MAF to coding variants in bona fide lymphoma driver genes before creating the report
- Driver gene list and coding variant classes are configurable (same options as maf_diff)
- Bundles the LLMPP curated tier-1 gene list as the default driver gene reference
- Produces a self-contained HTML report per tumour/normal pair with tumour and normal BAM tracks
