# Changelog

All notable changes to the `mfr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-11

### Fixed

- `envs/mfr.yaml` pinned `python=3.9`, which under `--use-conda` caused `_mfr_extract_chrom` jobs to
  fail with `ImportError: Unable to import required dependencies: numpy` before
  `extract_chrom_maf.py` (which never imports pandas/numpy itself) reached its own first line --
  visible as an empty `.extract.log` and a traceback pointing into the *driver* env's
  (e.g. GAMBL's `gambl-opv12`) site-packages. Snakemake's `script:` execution makes its own
  `snakemake` package importable inside an isolated conda env by grafting the driver env's
  site-packages onto that env's `sys.path` for the auto-generated preamble that runs before the
  user script; `gambl-opv12` is Python 3.8, so a Python 3.9 interpreter loading 3.8-compiled numpy
  C extensions via that graft failed on the ABI mismatch. Pinned to `python=3.8` to match.

### Changed

- `cluster_foci.R` now splits each chromosome's unique positions into gap-delimited chunks
  (wherever two consecutive positions are more than `h_max` apart) and clusters each chunk
  independently, instead of running one `dist()`/`hclust()` over every unique position on the
  chromosome. This is exact (not approximate) for any cut height `<= h_max`, since points on
  opposite sides of such a gap can never be merged into the same cluster below that height under
  the supported linkage methods -- but it bounds the O(n^2) `dist()` blowup to the largest chunk
  rather than the whole chromosome, which matters once a sample_set pools ~700 WGS samples onto
  one chromosome. Not a hard memory ceiling: a region with no gap `> h_max` anywhere (e.g. a
  shared hotspot) still forms a single large chunk.
- Renamed from the `mutation_foci` prototype to `mfr`; version reset to 1.0 as a fresh identity.
- Replaced the single genome-wide `inputs.master_maf` (all samples) with per-sample,
  bgzip + tabix-indexed `inputs.sample_maf`. A sample_set's jobs now only ever open the MAFs of
  samples actually in that sample_set.
- Replaced `prepare_maf.R` (full master-MAF read per sample_set) and the full-file
  read-then-filter-to-one-chromosome step in `cluster_foci.R` with a single new
  `_mfr_extract_chrom` rule (`src/python/extract_chrom_maf.py`) that streams each sample's MAF
  through `tabix {sample}.maf.gz {chrom}` one line at a time, dropping coding
  `Variant_Classification` rows inline, and writing straight through to the output file — never
  materializing a whole chromosome's rows as a table (no pandas/polars needed for this step).
- `cluster_foci.R` no longer filters by chromosome (its input is already scoped to one
  sample_set × chromosome by the upstream extraction rule) or reads a `chrom_column` option — it
  just does the clustering math on an already-small file.
- Removed `options.maf_sample_column` (was used to filter the master MAF by
  `Tumor_Sample_Barcode`; no longer needed since sample selection is now by file path, not by
  filtering a column inside a shared file).
- Single combined conda env / container image (`mfr`: R clustering deps + Python + htslib) used
  by both the extraction and clustering rules — nothing requires them to run in separate
  environments, so one image keeps the module and its lcr-scripts container simpler than
  maintaining two.

### Added

- Initial implementation of `mfr` module (as `mutation_foci`)
- Hierarchical clustering of non-coding mutation positions into "foci", scattered by sample_set
  and chromosome, selecting the cut height that maximizes mean silhouette width
