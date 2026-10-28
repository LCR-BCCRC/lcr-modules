# Changelog

All notable changes to the `vquest` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1] - 2026-07-15

This update was authored by Laura Hilton.

- **Per-request HTTP timeout:** `run_vquest.py` now patches `requests.post` with
  a configurable timeout (default 120 s, set via `options.request_timeout`) before
  calling the `vquest` library. Previously, a stalled TCP connection to
  `www.imgt.org` would hang the job indefinitely because the library makes no
  timeout provision of its own.

- **Concurrency limiting:** Running many `_vquest_run` jobs simultaneously causes
  IMGT V-QUEST to refuse connections (errno 111) at the TCP level. Each job
  declares `vquest: 1` in its `resources` block; pass `--resources vquest=N`
  (recommended N = 5–10) to cap concurrent IMGT connections. This limit is
  enforced by Snakemake's scheduler and applies equally to local and cluster runs.
  The pipeline now **refuses to start** if `--resources vquest=N` is omitted or
  if N exceeds `options.max_imgt_connections` (default: 10).

## [1.0] - 2026-05-25

This release was authored by Laura Hilton.

IMGT V-QUEST annotation module. Submits per-chain FASTA files to the IMGT V-QUEST
web service via the `vquest` Python package and writes AIRR-format TSV output, with
optional N-linked glycosylation annotation and merging back into a user-supplied
source TSV.

- Accepts per-chain FASTA from `igseqr` 1.0 or `mixcr` 1.3 (same wildcard
  convention as `igblast` 1.0: `{seq_type}`, `{sample_id}`, `{chain}`).
- Supported chains: IGH, IGK, IGL, IGKL (IG) and TRA, TRB, TRD, TRG (TR).
- Sequences longer than 10,000 bp are silently dropped before submission, as
  V-QUEST cannot process them.
- Species and molecule type are configurable via `options.species` and
  `options.molecule_type`; the `vquest` package batches sequences automatically
  (50 per request) so arbitrarily large inputs are handled transparently.
- Configurable retries (`options.vquest_retries`) absorb transient IMGT server
  errors without failing the job.
- Output: `{sample_id}.{chain}.vquest_airr.tsv` in AIRR Community format.
- Annotates N-linked glycosylation sites (`_vquest_annotate_glycosylation`): uses
  ANARCI to assign IMGT unique numbering and identifies NxS/T motifs (x ≠ Pro)
  that are present in the query but absent from the germline (SHM-acquired sites).
- Merges AIRR annotation and glycosylation results into a user-supplied source TSV
  (`_vquest_merge_final`), keyed on `options.source_id_column`; symlinked to
  `99-outputs/tsv/` as `{sample_id}.{chain}.merged.tsv`.
- Empty input FASTA and TSV files are handled gracefully throughout the pipeline.
- Conda env YAMLs live in the shared `envs/` directory at the repo root.

**Requirements:**

- Outbound HTTPS access to `www.imgt.org` from the executing node.
  On clusters where compute nodes are firewalled, run this module on a
  login or submit node, or use a workflow proxy.
- conda env: `vquest` installed from source via post-deploy script;
  no pre-built biocontainer is available — use `--use-conda`.
