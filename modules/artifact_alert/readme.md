# 🚨 Artifact Alert — Module README 🧬

A Snakemake workflow that builds a position-specific background mutation-rate index, ie how many times a non-reference aberration appears at a genomic position. Assists in alerting you to positions with recurrent technical artifacts (and real mutations) in NGS data. Can be run on capture panels small and large as well as WES and WGS bam files, with or without a .bed file (but .bed file recomended for faster execution). 

## 📋 Overview
Artifact Alert computes per-sample, per-position mutation rates (REF -> A, T, C, G, INS, DEL) and aggregates them across a cohort to produce a bgzip-compressed, tabix-indexed background mutation-rate index. The index is suitable for fast genomic lookups when filtering likely technical artifacts from variant calls.

Key points:
- Per-sample mutation-rate tables are computed, then aggregated across samples.
- Supports incremental updates: you can add new samples into an existing index.
- Per-chromosome parallel processing for speed; each chromosome is position-sorted before concatenation.
- Processed samples are tracked in a `sample_tracker` file to avoid duplication.
- can trigger remaking index from scratch by deleting the sample tracking file, or setting `reset_mutation_index` to `True` in the config.
- The final output index file is given to snakemake as a param NOT as an output--so snakemake will never delete the index if the pipeline is terminated prematurely or if it wants to remake the index.
- if a bed file is not provided, mpileup will run with no target regions.

## 🔬 Workflow Steps

### 1️⃣ Generate Pileup (`generate_pileup`)
Creates a pileup for each sample at targeted positions.
- Uses `samtools mpileup`.
- Filters by mapping and base quality.
- Runs only on regions in .bed file if one provided, else runs on whole bams.
- Output per-sample: `{sample_id}.pileup` (under `01-pileup/`).

### 2️⃣ Calculate Mutation Rates (`calculate_mutation_rates.py`)
Computes position-specific mutation rates for each sample.
- Calculates per-position aberration rates, aka how many times a non-reference base appears: `A`, `T`, `C`, `G`, `INS`, `DEL`.
- Applies minimum depth and other filters.
- Outputs per-sample mutation rate table: `{sample_id}_mutation_rates.tsv[.gz]` (under `02-mutation_rates/`).

### 3️⃣ Aggregate Mutation Rates (`aggregate_mutation_rates.py`)
Combines per-sample tables to create the cohort-wide background profile.
- Aggregates mean, standard deviation, and `M2` (see below) for each mutation type at every position.
- Supports reading an existing bgzip+tabix aggregated file and merging new samples without reprocessing old ones.
- Writes one bgzip-compressed, tabix-indexed file: `background_mutation_rates.tsv.gz` (under `03-aggregated/`).

## 📂 Output structure
- `artifact_alert/1.0/{PANEL_NAME}/`
  - `01-pileup/`
    - `{sample_id}.pileup`
  - `02-mutation_rates/`
    - `{sample_id}_mutation_rates.tsv` or `.tsv.gz`
  - `03-aggregated/`
    - `background_mutation_rates.tsv.gz` (bgzip + tabix) — final index
  - `logs/`
  - `{PANEL_NAME}_MutRateIndex_Sampletracker.txt` (tracks processed sample IDs)

## Usage (minimal example)
Provide a sample table with at least a `sample_id` column and include the Snakemake module:

```/dev/null/usage_example.py#L1-9
configfile: "path/to/config.yaml"

# Provide a table with at least a `sample_id` column
sample_table = pd.read_csv("samples.tsv")
config["lcr-modules"]["_shared"]["samples"] = sample_table

include: "artifact_alert/1.0/FetchMutationRate.smk"

rule execute_mut_rate_pipe:
    input:
        str(rules.all_artifact_alert.input)
```

## Output format

Columns present in `background_mutation_rates.tsv.gz`
Column | Description
---|---
`chromosome` | Chromosome name
`position` | Genomic position (1-based)
`ref_base` | Reference base
`n_samples` | Number of samples contributing to that position
`A_mean`, `A_std`, `A_M2` | Mean, std dev, and M2 for A errors
`T_mean`, `T_std`, `T_M2` | Mean, std dev, and M2 for T errors
`C_mean`, `C_std`, `C_M2` | Mean, std dev, and M2 for C errors
`G_mean`, `G_std`, `G_M2` | Mean, std dev, and M2 for G errors
`INS_mean`, `INS_std`, `INS_M2` | Mean, std dev, and M2 for insertion errors
`DEL_mean`, `DEL_std`, `DEL_M2` | Mean, std dev, and M2 for deletion errors

## Why we include `*_M2`
- `M2` is the running sum of squared deviations from the mean used by Welford’s algorithm.
- From `M2` you compute the sample variance as: variance = M2 / (n - 1) (for n ≥ 2) and `std = sqrt(variance)`.
- Storing `M2` makes incremental updates exact and numerically stable: you can merge an existing aggregated record with new samples without re-reading all original raw data.

## How on-the-fly aggregation is done (Welford)
Per-position objects are instantiated in `aggregate_mutation_rates.py` that track the running mean and variation for each mutation type. When you add a single new sample value `x` for a position, Welford’s online update is applied so mean and `M2` are updated without storing all individual sample values.
