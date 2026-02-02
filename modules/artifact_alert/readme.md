# 🚨 Artifact Alert 🧬

A Snakemake workflow for calculating the background mutation rates across capture panels to identify genomic positions prone to technical artifacts.

## 📋 Overview

Artifact Alert generates position-specific background error profiles from a cohort of samples. By calculating mutation rates at each genomic position across multiple samples, it helps distinguish true somatic variants from recurrent technical artifacts in NGS data.

This workflow creates a **background mutation rate index** that captures these systematic errors, allowing you to filter out likely artifacts in your variant calls.

## 🔬 Workflow Steps


### 1️⃣ Generate Pileup (`generate_pileup`)
Creates a pileup file for each sample at targeted genomic positions.
- Uses `samtools mpileup`
- Filters by mapping quality and base quality
- Restricts analysis to regions defined in target BED file

### 2️⃣ Calculate Mutation Rates (`calculate_mutation_rates`)
Computes position-specific mutation rates for each sample.
- Calculates mean error rate for each alternate base (A, T, C, G, INS, DEL)
- Applies minimum depth filtering
- Outputs per-sample mutation rate tables

### 3️⃣ Aggregate Mutation Rates (`aggregate_mutation_rates`)
Combines all samples to create a cohort-wide background mutation profile.
- Calculates mean and standard deviation across samples
- Produces final **background mutation rate index**
- One row per genomic position with summary statistics

## 📂 Output Structure

- `artifact_alert/1.0/{PANEL_NAME}/`
  - `01-pileup/`
    - `{sample_id}.pileup`
  - `02-mutation_rates/`
    - `{sample_id}_mutation_rates.tsv`
  - `03-aggregated/`
    - `background_mutation_rates.tsv` ⭐ **(Final output!)**
  - `logs/`

# Usage

The only required column for this module is "sample_id"

```
import pandas as pd

configfile: path/to/config

sample_table = pd.read_csv(input_file)

config["lcr-modules"]["_shared"]["samples"] = sample_table

include: ".../artifact_alert/1.0/FetchMutationRate.smk"

rule execute_mut_rate_pipe:
    input:
        str(rules.FetchMutationRate.input)
```

# Output format

Column | Description |
|--------|-------------|
| `chromosome` | Chromosome name |
| `position` | Genomic position |
| `ref_base` | Reference base |
| `n_samples` | Number of samples with coverage |
| `A_mean`, `A_std` | Mean and std dev for A errors |
| `T_mean`, `T_std` | Mean and std dev for T errors |
| `C_mean`, `C_std` | Mean and std dev for C errors |
| `G_mean`, `G_std` | Mean and std dev for G errors |
| `INS_mean`, `INS_std` | Mean and std dev for insertion errors |
| `DEL_mean`, `DEL_std` | Mean and std dev for deletion errors
