# Changelog

All notable changes to the `igv_reports` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.0] - 2026-07-09

### Added

- Initial implementation of `igv_reports` module wrapping [igvteam/igv-reports](https://github.com/igvteam/igv-reports) v1.16.0
- Accepts an arbitrary number of MAF inputs via `inputs.mafs` (a dict of label → path pattern); each label becomes a `{tool}` wildcard so one module include handles multiple callers
- Filters each input MAF to coding variants in bona fide lymphoma driver genes before creating the report, keeping reports focused and file sizes manageable
- Bundled LLMPP curated tier-1 gene list used by default (same as `maf_diff`); driver gene list and coding variant classes are fully configurable
- Uses `--fasta` with the `reference_files` genome FASTA rather than `--genome` to avoid chromosome naming mismatches between grch37 (no `chr` prefix) and UCSC-style genome identifiers
- Gencode annotation GTF added as an explicit track so gene annotations are embedded in the report and visible without network access
- Tumour and normal BAM index symlinks created with both `.bai` and `.crai` extensions pointing to the same source index, allowing pysam to discover whichever extension it expects regardless of alignment file type
- Self-contained HTML output (one file per tumour/normal pair per tool label) suitable for sharing and offline viewing
