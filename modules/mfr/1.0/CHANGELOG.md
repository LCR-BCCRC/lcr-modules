# Changelog

All notable changes to the `mfr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added

- `options.clustering_method`'s value is now baked into every `_mfr_cluster`/`_mfr_aggregate`/
  `_mfr_output_tsv` output path as a literal directory segment (e.g.
  `03-foci/dynamic/{sample_set}/...`), not a Snakemake wildcard (it's a single fixed value for the
  whole run, not scattered per-job) -- so switching methods and re-running no longer silently
  overwrites the previous method's results. `_mfr_input_maf`/`_mfr_input_subsets`/
  `_mfr_extract_chrom` are unaffected and stay shared across methods, since chromosome extraction
  doesn't depend on which clustering method is configured.
- `cluster_foci.R` now supports three interchangeable clustering methods, selected via the new
  `options.clustering_method` option (`"silhouette"`, `"density"`, or `"dynamic"`):
  - `"silhouette"`: the original mean-silhouette-width criterion.
  - `"density"`: a density/fragmentation-penalized score (one height per gap-chunk, like
    silhouette, but rewards tightly-packed clusters and penalizes excessive fragmentation).
  - `"dynamic"`: `dynamicTreeCut::cutreeDynamic()` (method = `"hybrid"`), which resolves cluster
    boundaries per dendrogram branch instead of picking one height for the whole chunk, so
    different foci in the same gap-chunk can end up with different effective cut heights.
- `options.gap` is now a standalone chunk-splitting parameter, decoupled from the height ceiling
  each method actually searches (`options.h_max` for silhouette/density, `options.cut_height` for
  dynamic). Previously `h_max` served double duty as both the chunk boundary and the search
  ceiling, so raising it for memory reasons also silently loosened what counted as a biologically
  valid focus. A runtime warning fires if `gap` is set below the active method's own ceiling,
  since gap-based chunking's exactness proof requires `gap >= h_max` (silhouette/density); this
  proof does not extend to the dynamic method, for which chunking is a pragmatic memory bound
  rather than a proven-exact optimization (see `cluster_foci.R`'s header for the detailed proof
  and caveats, including the `ward.D2` linkage exception).
- `options.min_cluster_size`, `options.deep_split`: `cutreeDynamic()` tuning (dynamic method only).
- `options.max_abs_core_scatter`, `options.min_abs_gap` (dynamic method only, optional): bound
  cluster tightness in actual bp rather than by point count, so an isolated pair of mutations far
  apart with nothing bridging them can be rejected even though `min_cluster_size` alone would
  accept any pair regardless of distance.
- `options.make_plots`: set `False` to skip diagnostic plot generation (a cheap placeholder is
  still written so declared outputs exist).
- Second diagnostic plot output, `{chrom}.h_score.pdf` (silhouette/density methods only): one
  point per gap-chunk, its chosen height against the score achieved there. The dynamic method has
  no native scalar score to plot (a placeholder is written instead) — retroactively scoring its
  output with e.g. silhouette width would unfairly favour whichever method is actually optimizing
  that metric.
- `schemas/base-1.0.yaml`: symlinked to the shared `schemas/base/base-1.0.yaml`, matching
  repo-wide convention (was previously missing for this module).

### Changed

- `{chrom}.silhouette.pdf` renamed to `{chrom}.h_distribution.pdf`, since it's no longer
  silhouette-specific and its content changed from a per-chunk score-vs-height curve (unreadable
  once a chromosome has more than a handful of gap-chunks) to a per-chromosome histogram of chosen
  cut heights (silhouette/density: one point per gap-chunk; dynamic: one point per focus, since
  that method has no single height per chunk).
- Cluster ids within a gap-chunk are now renumbered by ascending minimum position after
  clustering. `cutree()`'s (and especially `cutreeDynamic()`'s) raw label order is not guaranteed
  to follow the left-to-right spatial order of the input positions, which could previously produce
  non-monotonic `group` values within a chunk (e.g. row 5 in a higher-numbered group than row 6).
- Default `hclust_method` changed from `"centroid"` to `"average"`: centroid/median linkage can
  produce non-monotonic merge heights, which is a worse fit for `cutreeDynamic`'s branch-jump
  heuristics (which assume roughly monotonic merge heights) than it is for the single-height
  `cutree()` methods.
- Histogram bars in diagnostic plots no longer have a white outline (`colour = "white"` ->
  `colour = NA`): with `binwidth = 1` spread across a wide height range, bars can be a fraction of
  a point wide, and the white outline was rendering as a blank white plot instead of visible bars.

## [1.0] - 2026-07-11

### Changed

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
