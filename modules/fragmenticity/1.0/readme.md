# fragmenticity 

## Overview
Fragmenticity produces fragmentation‑score (FS) outputs from cfDNA BAMs via two Snakemake workflows:
- regional scoring (`fragmenticity.smk`) — per‑sample genome‑wide + BED‑regional FS and histogram.
- FS atlas creation (`ProfileFragments.smk`) — builds the fragmentation‑score reference (length 1–900) from tumor vs healthy profiles.

## Included scripts
- `fragmenticity.smk` — Snakemake ruleset to run regional scoring across samples (rules: `Generate_FS_Regional`, `Generate_FS_all`).
- `ScoreFragsRegional.py` — regional FS engine (genome‑wide sampling, per‑region FS, adaptive FS, PNG histogram).
- `ProfileFragments.smk` — workflow to aggregate profiles and compute FS atlas (final rule: `all_CreateFS`).
- `ProfileMutatedReads.py` — extract mutated read length profiles from BAM + MAF.
- `ProfileHealthyReads.py` — sample healthy read‑length profiles from BAMs.
- `CalculateFSDistribution.py` — bootstrap + compute log2 FS per length (1–900).

## Inputs & assumptions
- BAMs: scripts expect UMIs under tag `MI` (reads missing `MI` are skipped).
- BED for regional scoring: 4 columns (`chr`, `start`, `end`, `name`); regions with the same `name` are collapsed.
- Fragment score reference (used by `ScoreFragsRegional.py`): either a single column (no header) listing scores for lengths 1..N OR two columns with headers `length` and `fragmentation_score`.
- Default notable parameters (can be set via config or CLI in scripts): `read_count` ≈ 100000, `min_reads_per_region` ≈ 1000, `top_k_regions` ≈ 10, length range 50–900bp.

## Outputs (per sample)
- `{sample}_FragmentScore.tsv` — main summary (genome‑wide FS, adaptive FS, max regional FS, mean length, counts).
- `{sample}_RegionalFragmentScores.tsv` — per‑region FS table (if regions had sufficient reads).
- `{sample}_regional_fragment_distribution.png` — fragment length distribution plot.
- FS atlas: `FragmentScoreAtlas/99-FSDistribution/FS_distribution.tsv` (per `ProfileFragments.smk`).
