# Changelog

All notable changes to this module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-06-29

### Added
- Initial release of the WGCNA co-expression module 
- `_wgcna_normalize_expression`: DESeq2 VST normalization and limma batch correction (raw mode)
- `_wgcna_filter_variance_genes`: per-gene median and MAD variance filtering
- `_wgcna_get_coexpression_modules`: soft-threshold power selection and `blockwiseModules` WGCNA run
- Support for both `raw` (count matrix) and `normalized` (pre-normalized matrix) input modes
- Per-pathology branching; set `pathologies: null` to run over all samples as a single cohort
- Container support via Docker/Apptainer (`ghcr.io/LCR-BCCRC/lcr-scripts/wgcna:1.73`)
- Diagnostic plots written to `{pathology}/plots/` (PCA, distribution histograms, dendrogram)