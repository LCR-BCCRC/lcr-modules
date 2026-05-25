# Changelog

All notable changes to the `vquest` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-05-25

This release was authored by Laura Hilton.

Initial implementation of the IMGT V-QUEST module.

Submits per-chain FASTA files to the IMGT V-QUEST web service via the
`vquest` Python package (v0.0.10) and writes AIRR-format TSV output.

- Accepts per-chain FASTA from igseqr 1.0 or mixcr 1.3 (same wildcard
  convention as igblast 1.0: `{seq_type}`, `{sample_id}`, `{chain}`)
- Supported chains: IGH, IGK, IGL, IGKL (IG) and TRA, TRB, TRD, TRG (TR)
- Species is configurable via `options.species` (default: `human`)
- The `vquest` package batches sequences automatically (50 per request)
  so arbitrarily large FASTA inputs are handled transparently
- Output: `{sample_id}.{chain}.vquest_airr.tsv` in AIRR Community format

**Requirements:**

- Outbound HTTPS access to `www.imgt.org` from the executing node.
  On clusters where compute nodes are firewalled, run this module on a
  login or submit node, or use a workflow proxy.
- conda env: Python ≥3.6, `vquest==0.0.10` (installed via pip);
  no pre-built biocontainer is available — use `--use-conda`.
