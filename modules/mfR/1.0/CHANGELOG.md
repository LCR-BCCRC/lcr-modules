# Changelog

All notable changes to the `mfR` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-11

### Changed

- Renamed from the `mutation_foci` prototype to `mfR`; version reset to 1.0 as a fresh identity.
- Replaced the single genome-wide `inputs.master_maf` (all samples) with per-sample,
  bgzip + tabix-indexed `inputs.sample_maf`. A sample_set's jobs now only ever open the MAFs of
  samples actually in that sample_set.
- Replaced `prepare_maf.R` (full master-MAF read per sample_set) and the full-file
  read-then-filter-to-one-chromosome step in `cluster_foci.R` with a single new
  `_mfR_extract_chrom` rule (`src/python/extract_chrom_maf.py`) that streams each sample's MAF
  through `tabix {sample}.maf.gz {chrom}` one line at a time, dropping coding
  `Variant_Classification` rows inline, and writing straight through to the output file — never
  materializing a whole chromosome's rows as a table (no pandas/polars needed for this step).
- `cluster_foci.R` no longer filters by chromosome (its input is already scoped to one
  sample_set × chromosome by the upstream extraction rule) or reads a `chrom_column` option — it
  just does the clustering math on an already-small file.
- Removed `options.maf_sample_column` (was used to filter the master MAF by
  `Tumor_Sample_Barcode`; no longer needed since sample selection is now by file path, not by
  filtering a column inside a shared file).
- Split the module's single conda env / container image into two: `r_foci` (clustering) and
  `python_tabix` (extraction; Python + htslib only, no pandas). Correspondingly, the lcr-scripts
  side now builds two images (`mfr_r`, `mfr_tabix`) instead of one `mutation_foci` image.

### Added

- Initial implementation of `mfR` module (as `mutation_foci`)
- Hierarchical clustering of non-coding mutation positions into "foci", scattered by sample_set
  and chromosome, selecting the cut height that maximizes mean silhouette width
