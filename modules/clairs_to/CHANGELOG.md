# Changelog

All notable changes to the `clairs_to` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2026-09-14

This release was authored by Giuliano Banco.

Variant calling is now parallelized by chromosome instead of one whole-genome `run_clairs_to` invocation per sample. A checkpoint reads the reference chromosome list, `_clairs_to_call_variants` runs once per chromosome via ClairS-TO's `--bed_fn`, and the resulting per-chromosome VCFs are merged (`bcftools concat`/`sort`/index) back into the same final path the old single-shot rule used to produce, so `_clairs_to_filter` and everything downstream needed no changes.

`--bed_fn` is ClairS-TO's only region-restriction mechanism, so it's reused for two purposes: restricting to a target panel (`options.target_regions_bed` — leave `""` for unrestricted whole-chromosome calling), and the chromosome split itself. Every per-chromosome job gets its own BED, either the target BED intersected with that chromosome, or (if none given) that chromosome's own full extent read from `genome.fa.fai`. Do not also pass `--bed_fn` via `clairs_to_args` — it will conflict with the one the module adds automatically.

The chromosome list itself defaults to `reference_files`' `main_chromosomes.txt` (chr1-22 + chrX, no chrY), overridable via `options.chromosomes_file` if you need chrY or a different custom contig set called.

`run_clairs_to`'s own `--output_dir` (raw per-chromosome scratch: `tmp/`, `tmp_TUMOR/`, `logs/`, the pre-combine SNV/indel VCFs) can be redirected onto different storage via `options.intermediate_results_dir_base`, e.g. to avoid heavy small-file traffic on a disk that shouldn't take it. This must be shared/network-visible storage, not node-local (`/tmp`, `/var/tmp`) — per-chromosome jobs and their cleanup jobs can land on different compute nodes, and a node-local path is invisible from any node other than the one that wrote it, which would silently break cleanup. The final combined per-chromosome VCF always stays in the module's own managed directory regardless of this setting, so downstream rules are unaffected by where scratch space is redirected.

`threads`/`resources` config keys changed: `clairs_to` (previously sized for one whole-genome job per sample) is now `clairs_to_run` (sized per chromosome — expect ~23 concurrent jobs per sample for hg38's default contig list), and a new `clairs_to_merge_vcfs` key was added for the merge step.

## [1.0] - 2026-02-20

This release was authored by Giuliano Banco.

This module is designed to only work with a sample table that includes the columns 'chemistry' and 'platform'. These are required by schemas. Note that ClairS-TO was developed to support R10 data but not R9, however, you can specify R10 as the chemistry value and use the R10 models.

This module only works with PromethION and hg38 data, as of version 1.0.