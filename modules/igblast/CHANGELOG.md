# Changelog

All notable changes to the `igblast` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-05-20

This release was authored by Laura Hilton.

Initial implementation of IgBLASTn (v1.17.1) as a standalone module for aligning
receptor sequences to IMGT germline databases.

- Accepts per-chain FASTA input from either `mixcr` 1.3 or `igseqr` 1.0 via a
  user-configured `sample_fasta` path.
- Supported chains: IGH, IGK, IGL, TRA, TRB, TRD, TRG (standard); IGKL for
  combined kappa/lambda FASTA output from igseqr.
- Runs `igblastn` with IMGT germline V/D/J databases; outputs raw fmt7 as an
  intermediate and a human-readable TSV as the final deliverable.
- TSV columns: sequence_id, productive, top_v_allele, top_v_identity,
  mutated_status, top_v_bit_score, top_v_evalue, top_j_allele, top_v_alleles,
  all_v_alleles, all_v_identities, all_mutated_status, all_bit_scores,
  all_evalues, btop_mutations.
- igBLAST post-processing specific to MiXCR (merging results back into clonotype
  tables) is handled by `mixcr` 1.3 via its `run_igblast_postprocessing` option.
