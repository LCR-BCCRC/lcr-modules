# Changelog

All notable changes to the `maf_diff` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-09

### Added

- Initial implementation of `maf_diff` module
- Compares MAF files from two variant callers per tumour/normal pair
- Outputs caller-specific MAFs for variants unique to each caller (all original columns preserved)
- Matching is performed on configurable key columns (default: Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
- Long-format summary TSV with one row per category (intersect, caller1-only, caller2-only) and columns for variant counts, driver variant counts, and driver gene lists
- Bundled LLMPP curated tier-1 lymphoma driver gene list (`etc/any_tier1_BL_FL_DLBCL_MCL_MZL_PMBL.tsv`) used by default for driver gene statistics
- Driver gene analysis restricted to coding variant classes (missense, nonsense, frameshift, in-frame indels, splice site, translation start, nonstop); configurable via `options.coding_variant_classes`
- Driver gene reporting can be disabled by setting `options.driver_genes: ""`
