# Changelog

All notable changes to the `deepsomatic` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-09-14

This release was authored by Giuliano Banco.

This module calls somatic variants from ONT long-read data using Google's DeepSomatic. It supports two calling modes via `options.calling_mode`: `unmatched` (tumour + an unmatched normal, the default) and `tumor_only`. The sample table must include a `chemistry` column (`R9` or `R10`); `options.normal_name` maps chemistry to the correct unmatched-normal sample per run. As of writing, this module only works with hg38 data.

Variant calling is parallelized by chromosome: a checkpoint reads the reference chromosome list (`reference_files`' `main_chromosomes.txt`, chr1-22 + chrX, no chrY by default — override via `options.chromosomes_file` if you need chrY or a different contig set), `_deepsomatic_call_variants` runs once per chromosome via DeepSomatic's `--regions`, and the resulting per-chromosome VCFs are merged (`bcftools concat`/`sort`/index) into one final VCF before filtering/gnomAD annotation.

`options.target_regions_bed` optionally restricts calling to a target panel — leave `""` for unrestricted whole-chromosome calling. It's a dedicated option rather than something to bake into `deepsomatic_args`: the module already appends its own `--regions` per chromosome, and DeepSomatic treats multiple `--regions` values as a union rather than an intersection, so simply concatenating a chromosome name with a whole-genome BED wouldn't restrict a chromosome's job to just its own targeted subregions. When a target BED is set, it's intersected with each chromosome automatically.

DeepSomatic always writes intermediate/scratch files somewhere — `--intermediate_results_dir` is passed on every run and defaults to the module's own managed `deepsomatic_temp/` subdirectory. Each chromosome's intermediate directory is deleted (`_cleanup_intermediate_dir`) as soon as that chromosome's own VCF is written, not after the whole sample finishes, unless `options.cleanup_toggle` is set to `false`. The base directory for `--intermediate_results_dir` can be redirected via `options.intermediate_results_dir_base` (e.g. to avoid heavy small-file traffic on a disk that shouldn't take it) — this must be shared/network-visible storage, not node-local (`/tmp`, `/var/tmp`): per-chromosome jobs and their cleanup jobs can land on different compute nodes, and a node-local path is invisible from any node other than the one that wrote it, which would silently break cleanup.
