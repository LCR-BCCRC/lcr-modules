# Changelog

All notable changes to the `panel_of_normals` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- added an optional GATK section (`options: gatk: enabled`, default False) that builds the read-count panel of normals (CreateReadCountPanelOfNormals), padded target intervals, GC annotation and common gnomAD SNP sites for the `gatk_cnv` module, under `99-outputs/gatk/{seq_type}--{genome_build}/{capture_space}/`
- the GATK panel for each combo uses the same samples as the CNVkit and PureCN panels (the module's samples table, without `flat_ref_combos`), so the curated selection, including any tumours chosen as pseudo-normals, is reused
- disabled by default, and the new subdirectories are numbered after the existing ones, so existing outputs and job plans are unchanged (checked by dry run: identical jobs and output paths with the section disabled)

## [1.0] - 2026-02-20

This release was authored by sgillis.

Initial set up of module. It contains the relevant parts for generating panel of normals reference files for modules `cnvkit` and `purecn`.

Note that this module requires the `capture_space` wildcard and maintains it through all output file names.
