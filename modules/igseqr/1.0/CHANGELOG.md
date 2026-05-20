# Changelog

All notable changes to the `igseqr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-05-20

This release was authored by Ryan Morin.

Initial implementation of IgSeqR (v1.0.1) for de novo assembly of B cell tumor
immunoglobulin transcripts from RNA-seq data.

- Input: aligned BAM (mrna seq_type) + HISAT2 splice-aware index
- Output per sample:
  - `{sample_id}_IGH_transcripts.fasta` / `{sample_id}_IGKL_transcripts.fasta`
  - `{sample_id}_IGH_report.tsv` / `{sample_id}_IGKL_report.tsv`
  - `{sample_id}_IGH_dominant_report.tsv` / `{sample_id}_IGKL_dominant_report.tsv`
  - `{sample_id}_IGH_TPM_filtered.fasta` / `{sample_id}_IGKL_TPM_filtered.fasta`
- Install rule clones IgSeqR from GitHub and runs setup.sh once; sentinel file
  prevents redundant re-installation.
- Conda env provides: blast 2.13.0, hisat2 2.2.1, kallisto 0.48.0,
  samtools 1.16.1, trinity 2.13.2.
- Container support added (container_envs.igseqr points to ghcr.io image).
- Only mrna seq_type is supported; pairing_config runs tumours unpaired.
- New `_igseqr_get_hisat_ref` rule downloads the recommended GRCh38 HISAT2
  splice-aware index (grch38_snptran) automatically. The download URL is
  configurable via `hisat_ref_url` (a version-keyed dict in config); the index
  is stored under `{inputs_dir}/hisat_ref/{version}/` and the path is derived
  in the module — no manual `hisat_ref` config entry required.

Reference: Carrington et al. (2022) Nature Protocols. doi:10.1038/s41596-022-00700-6
