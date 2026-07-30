# Changelog

All notable changes to the `mhc_hammer` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-07-30

This release was authored by Ryan Morin.

- Ported the DNA analysis arm of MHC Hammer ([McGranahanLab/mhc-hammer](https://github.com/McGranahanLab/mhc-hammer)) to native Snakemake rules, rather than wrapping upstream's Nextflow pipeline. RNA allelic expression/imbalance/repression and alternative splicing are out of scope for this version (see README).
- Novoalign and HLA-HD are user-supplied filesystem paths (`options.novoalign_dir`, `options.hlahd_dir`), not conda/container-managed -- both are gated behind manual, registration/licence-restricted downloads that can't be automated into a build.
- MHC Hammer's own `bin/*.R`/`bin/*.sh` scripts are **not bundled in this module** and never will be: upstream is distributed under a Cancer Research Horizons academic-use licence that prohibits redistribution and modification (see README "Prerequisites"). `options.mhc_hammer_scripts_dir` points at a clone the user obtains and licenses themselves, following the same pattern already used by this repo's `lymphgen`/`dlbclass` modules for similarly-restricted third-party code.
- VEP is user-supplied (`options.vep_path`, `options.vep_cache`), following the `vcf2maf` module's own pattern, to avoid the bioconda/Perl `ensembl-vep` dependency conflicts already documented for that module.
- The MHC/IMGT reference bundle (versioned, non-gated Zenodo download) is fetched and cached by the module itself (`_mhc_hammer_download_reference`).
- Sample grouping uses standard `oncopipe` tumour/normal pairing (`CFG["paired_runs"]`) for WES pair-level rules, plus a custom `patient_id`-wildcarded pattern (mirroring `modules/gridss/2.0/gridss.smk`) for the two per-patient steps (HLA-HD typing, personalised reference construction) that are shared across all of a patient's tumour samples rather than recomputed per run.
- Requires exactly one germline WES sample per patient (errors otherwise) -- a v1 simplification of upstream's behaviour of tolerating multiple germline samples per patient by silently picking one.
- Rule shell commands calling `mhc_hammer_scripts_dir` scripts were cross-checked against a local read-only clone of upstream's `bin/*.R`/`bin/*.sh` (never copied into this repo -- read-only verification only, per the licensing note above). This caught and fixed several real bugs: upstream hardcodes the literal token `wxs` (not a generic seq_type) into filenames several cohort-level R scripts pattern-match on; `make_mutation_table.R` and `make_cohort_overview_table.R` expect to run from a flat directory of bare-basename symlinks to their inputs (mirroring Nextflow's per-task staging) and require a synthesised `inventory.csv`/`hlahd_germline_samples.txt` this module now generates (`_mhc_hammer_generate_inventory`, derived from `CFG["samples"]`/`CFG["paired_runs"]`/`purity_ploidy_file` since lcr-modules has no equivalent of upstream's own samplesheet); `hlahd_parse_output.R` writes into its `--hlahd_folder` argument, requiring an extra move step; `concatenate_dna_analysis_tables.R`'s `--genes` needed to be an unquoted (word-split) argument list, and its output filename didn't match what the module declared.
