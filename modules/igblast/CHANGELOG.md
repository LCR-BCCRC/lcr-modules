# Changelog

All notable changes to the `igblast` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-05-20

This release was authored by Laura Hilton.

Standalone IgBLASTn (v1.17.1) module for V(D)J annotation of IG/TCR sequences
against IMGT germline databases, with optional N-linked glycosylation annotation
and merging back into a user-supplied source TSV.

- Accepts per-chain FASTA input from `igseqr` 1.0 or `mixcr` 1.3 via a
  user-configured `sample_fasta` path.
- Supported chains: IGH, IGK, IGL, TRA, TRB, TRD, TRG (standard); IGKL for
  combined kappa/lambda FASTA from igseqr.
- Runs `igblastn` against IMGT germline V/D/J databases and writes AIRR Community
  TSV output (`{sample_id}.{chain}.igblast_airr.tsv`).
- Optional constant-region database support via `options.c_region_db`.
- Annotates N-linked glycosylation sites (`_igblast_annotate_glycosylation`): uses
  ANARCI to assign IMGT unique numbering and identifies NxS/T motifs (x ≠ Pro)
  that are present in the query but absent from the germline (SHM-acquired sites).
- Merges AIRR annotation and glycosylation results into a user-supplied source TSV
  (`_igblast_merge_final`), keyed on `options.source_id_column`; symlinked to
  `99-outputs/tsv/` as `{sample_id}.{chain}.merged.tsv`.
- Empty input FASTA and TSV files are handled gracefully throughout the pipeline.
- Conda env YAMLs live in the shared `envs/` directory at the repo root.
