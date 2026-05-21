# Changelog

All notable changes to the `igseqr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-05-20

This release was authored by Ryan Morin.

Initial implementation of IgSeqR (v1.0.1) for de novo assembly of B cell tumor
immunoglobulin transcripts from RNA-seq data.

The monolithic `igseqr` shell script has been decomposed into individual Snakemake
rules so that each step runs as a separate cluster job and compute-heavy steps can
execute in parallel across samples:

- `_igseqr_hisat2_align` — HISAT2 alignment, coordinate-sorted BAM + index
- `_igseqr_filter_reads` — extract unmapped reads and reads overlapping IGH/IGK/IGL
  loci; convert to paired FASTQ
- `_igseqr_trinity_assemble` — Trinity de novo assembly of filtered reads (once per
  sample, shared across chains)
- `_igseqr_blastn` — blastn of Trinity contigs against IMGT germline database (per
  chain); BLAST databases are bundled with IgSeqR and deployed by the post-deploy
  script to `$CONDA_PREFIX/bin/data/igseqr/IMGT/{chain}/`
- `_igseqr_extract_transcripts` — extract BLAST-hit contigs from Trinity FASTA with
  samtools faidx
- `_igseqr_kallisto_index` — build kallisto index from chain-specific transcripts
- `_igseqr_kallisto_quant` — kallisto quantification against IG-filtered reads
- `_igseqr_make_report` — generate `_report.tsv`, `_dominant_report.tsv`, and
  `_TPM_filtered.fasta` via `src/igseqr_report.py`

Config changes:

- `options.igseqr_run` replaced by `options.hisat2_align`, `options.trinity_assemble`,
  and `options.blastn` (extra flags for each tool)
- `scripts.igseqr_report` added pointing to the new report script
- `threads` and `resources` entries added for each new rule; `igseqr_run` entries
  removed

- Input: paired FASTQ (mrna seq_type) + HISAT2 splice-aware index
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
  splice-aware index automatically. The download URL is configurable via
  `hisat_ref_url` (a version-keyed dict in config); the index is stored under
  `{inputs_dir}/hisat_ref/{version}/` and the path is derived in the module —
  no manual `hisat_ref` config entry required.
- Inputs changed from BAM+BAI to a paired FASTQ (R1/R2); `sample_bam` and
  `sample_bai` config keys replaced by `sample_fastq_1` and `sample_fastq_2`.
- `{genome_build}` wildcard removed throughout; all output directories are now
  structured as `{seq_type}--{hisat_ref_version}` so outputs are labelled by
  the reference used rather than the genome build.

Reference: Carrington et al. (2022) Nature Protocols. doi:10.1038/s41596-022-00700-6
