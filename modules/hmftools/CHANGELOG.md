# Changelog

All notable changes to the `hmftools` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2021-12-29

- Updates to the GRIDSS-PURPLE-LINX pipeline, incorporating new and better ways of handling unmatched tumours. 
- Added two small, additive prerequisite rules for the new `neo` module (HMF Tools' NEO neoantigen prediction), deliberately kept separate from `_hmftools_linx`/`_hmftools_purple_matched` themselves rather than adding to those rules' own `output:` lists -- under `--rerun-triggers mtime`, adding a new required output to an *existing* rule would make every already-completed sample look incomplete and force a rerun (LINX, and worse, PURPLE-and-everything-downstream-of-it) just because a newly-declared file didn't exist yet for old runs:
    - `_hmftools_linx_neo_epitopes`: re-invokes LINX with the exact same real inputs `_hmftools_linx` already uses, plus `-write_neo_epitopes` (confirmed to already exist in the pinned `hmftools-linx=1.15` by downloading the closest real release jar, `linx-v1.16`, and finding the literal string compiled into `FusionDisruptionAnalyser.class` -- no LINX version bump needed), into its own output directory. Costs one redundant LINX invocation per pair that wants neo-epitope output, but leaves `_hmftools_linx`'s own already-completed output completely untouched for every other sample/user.
    - `_hmftools_output_purple_somatic_vcf`: PURPLE already writes a real `<tumour_id>.purple.somatic.vcf.gz` as a standard side-effect of being given `-somatic_vcf` (which `_hmftools_purple_matched` already does) -- it was simply never declared as a tracked output. This rule depends only on one of `_hmftools_purple_matched`'s own *already-tracked* outputs (`purity.tsv`) as a dependency anchor for correct ordering/staleness, then symlinks the adjacent, already-real file directly -- `_hmftools_purple_matched`'s own `output:` is completely unchanged.
    - Both verified via a real local dry-run against demo data, including a real simulation of "PURPLE and LINX already completed" (their own declared outputs pre-touched with correctly-ordered timestamps) confirming neither rule appears in the resulting job list -- only the two new rules and their own downstream consumers do.

## [1.0] - 2020-07-29

This release was authored by Laura Hilton.
- This module takes output from the `gridss` module to perform purity/ploidy/CNV calling (PURPLE) and generate SV plots (LINX). 
- Is able to run in unmatched normal mode, but the results will be substantially noisier than on matched samples. 
- For optimal pipeline efficiency, include both `modules/gridss/1.0/gridss.smk` and `modules/hmftools/1.0.smk` in the same Snakefile, with `_hmftools_all` as the target rule. 
