# drivemutr

# Purpose

The `drivemutr` is a Level 4 module that operates on a cohort-level `MAF` and copy-number matrix, together with the co-expression outputs of the `WGCNA` module, to use DriveMuTR and search **aSHM** (aberrant Somatic HyperMutation) and other user-supplied target regions for non-coding **regulatory driver** mutation foci — clusters of mutations whose presence is associated with a change in the expression of the host gene — for matched `genome, mrna` data. It generates `UCSC custom track` `TSV` files as outputs.

Unlike most modules, `drivemutr` is *cohort-level and per-gene* rather than per-sample: it consumes a single MAF + copy-number matrix spanning the whole cohort, scatters one job per target gene via a checkpoint, and gathers the per-gene results into UCSC custom tracks. Its working wildcards are therefore `{gene}` and `{lam}` (one hierarchical-grouping lambda value) rather than the usual `seq_type`/`sample_id`.

For each target gene the module:

1. Collects the somatic mutations (SSM/MAF) and copy-number values that fall in the region, and scores each mutation for deleteriousness with **CADD**.
2. Groups nearby mutations into **foci** by hierarchical clustering (the `lambda` parameter controls how aggressively positions are merged), optionally refined with a sliding window.
3. Links each mutated sample to its **WGCNA co-expression module** and host-gene expression, keeping only samples that have matched DNA *and* RNA.
4. Fits per-focus models — a module-eigengene residual linear model, a genomic model, and a **RuleFit + Shapley** interaction model — to find foci significantly associated with expression.
5. Runs **motifbreakR** TF-motif disruption analysis on the significant foci and tests, with a Fisher exact test, whether the disrupted motifs are enriched for transcription factors expressed in the pathology.

The per-gene results are then aggregated into UCSC custom track TSV files for genome-browser visualization.

## Output

All module outputs are written under `{root_output_dir}/drivemutr-1.0/`. The final UCSC custom tracks are symlinked into `99-outputs/ucsc_custom_tracks/`, with one TSV per track **per `lambda` value** (e.g. `lambda: [1]` produces the `_lambda_1` suffix; `lambda: [0.5, 1]` adds `_lambda_0_5`):

| File | Description |
|---|---|
| `Mutation_Points_{lam}.tsv` | Per-position mutation track. Colour encodes whether the sample is RNA-matched (blue) or DNA-only (red); score encodes the genomic-model p-value of the focus. |
| `Mutation_Blocks_{lam}.tsv` | Foci (clustered mutation blocks) as genomic ranges. |
| `Transcription_Factors_{lam}.tsv` | Disrupted TF binding motifs at the significant foci. |

Per-gene results and diagnostic plots are kept in the module's per-gene working directory (`03-per_gene/{gene}/`):

| File | Description |
|---|---|
| `final_results.rds` | Terminal per-gene result object (foci, models, TF results) — the input to the aggregation step. |
| `height_plot.pdf` | Hierarchical-clustering merge heights used to define foci. |
| `shap_plots.pdf` | RuleFit / Shapley interaction plots for the gene's foci. |
| `sanity_check_plot.pdf` | Boxplot of host-gene expression by mutation-foci group. |

Intermediate per-gene `.rds` files (one per pipeline step) are marked `temp()` and removed automatically once the gene's `final_results.rds` is produced.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    drivemutr:
        inputs:
            # Cohort-level single files (no wildcards)
            ssm_maf: "data/mutations.maf"                           # MAF (TSV) of simple somatic mutations
            cnv_matrix: "data/cn_matrix.tsv"                        # copy-number matrix (TSV): col 1 = sample_id, rest = per-gene CN
            wgcna_coexpression_modules: "results/WGCNA/coexpression_modules.rds"  # WGCNA module output
            wgcna_filtered_expression: "results/WGCNA/filtered_expression.tsv"    # WGCNA module output (samples x genes)
            tf_names_file: "data/TF_names.txt"                     # one TF gene symbol per line (e.g. Lambert et al. 2018)
        reference_params:
            genome_build: "grch37"                                 # grch37/hg19/37 or grch38/hg38/38
            pathology: "DLBCL"                                     # stored in the Pathology column of the output
            ashm_regions: True                                    # include the aSHM region set from GAMBLR.data
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"

configfile: "../modules/drivemutr/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/drivemutr/1.0/drivemutr.smk"

rule all:
    input:
        rules._drivemutr_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/drivemutr/CHANGELOG.md)