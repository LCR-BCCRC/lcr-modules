
# Changelog

All notable changes to the `deepsomatic` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-08-10

This release was authored by Giuliano Banco.

Initial release of the `deepsomatic` module, which calls somatic variants from Oxford
Nanopore Technologies (ONT) long-read data using DeepSomatic. The workflow calls variants,
indexes the VCF with tabix, filters against a panel of normals and quality/depth/VAF
thresholds, annotates with gnomAD population frequencies, and symlinks a final filtered
VCF into the outputs directory.

### Requirements and constraints
- Requires a sample table with `chemistry` and `platform` columns (enforced by schemas).
- Only compatible with PromethION and hg38 data as of version 1.0.
- Requires container usage for the DeepSomatic calling step. The bcftools steps (index,
    filter, gnomAD annotation) support either conda or container.

### Features
- Two calling modes, set via `calling_mode` in the config.
    - `unmatched` (uses an unmatched normal BAM with the ONT model)
    - `tumor_only` (no normal, ONT_TUMOR_ONLY model). Normal-based filters are only
        applied in unmatched mode.
- Maps tumour chemistry (R9 or R10) to a normal sample name via the `normal_name`
    config option.
- Filters variants to gnomAD AF < 0.0001 (missing AF is treated as 0).
- Optional cleanup of DeepSomatic intermediate files with `cleanup_toggle` (recommended).