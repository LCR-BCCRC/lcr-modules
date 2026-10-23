# WGCNA Co-expression Module

## Purpose

The `WGCNA` module is a Level 3 module that operates on gene expression data to identify co-expressed gene modules using [WGCNA](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559) (Weighted Gene Co-expression Network Analysis). It supports two input modes:

| Mode | Input | Steps |
|---|---|---|
| `raw` | Raw count matrix (genes × samples, TSV) | DESeq2 VST normalization → limma batch correction → variance filtering → WGCNA |
| `normalized` | Pre-normalized matrix (genes × samples, RDS or TSV) | Variance filtering → WGCNA |

**Final outputs** (symlinked into `results/wgcna/1.0/02-outputs/{pathology}/`):

| File | Description |
|---|---|
| `coexpression_modules.rds` | Named vector mapping each gene to its module colour |
| `coexpression_modules.tsv` | The same gene → module assignments as a flat table (`gene`, `module` columns) |
| `network.rds` | Full WGCNA network object from `blockwiseModules` |

**Intermediate files** (written to `results/wgcna/1.0/01-wgcna/{pathology}/`):

| File | Description |
|---|---|
| `normalized_expression.tsv` | Batch-corrected VST matrix (genes × samples, with an `hgnc_symbol` column); raw mode only |
| `filtered_expression.tsv` | Variance-filtered expression matrix (samples × genes) passed to WGCNA |

> **Normalized mode note:** when `input_mode: "normalized"` is used there is no per-pathology split. All outputs land under the fixed path `results/wgcna/1.0/02-outputs/user_provided/`.

> **Raw mode, no pathologies:** set `pathologies: null` (or leave the list empty) to normalize over all samples without pathology filtering — outputs then land under `results/wgcna/1.0/02-outputs/all/`.


**Diagnostic plots** (written to `results/wgcna/1.0/01-wgcna/{pathology}/plots/`):

| File | Description |
|---|---|
| `PCA_uncorrected_{variable}.pdf` | PCA coloured by each batch/bio variable before correction (raw mode) |
| `PCA_corrected_{variable}.pdf` | PCA coloured by each batch/bio variable after correction (raw mode) |
| `Primary_Distribution.pdf` | Histogram of all expression values before filtering |
| `Median_Distribution_Before_Filtering.pdf` | Per-gene median distribution before filtering |
| `MAD_Distribution_Before_MAD_After_Median_Filtering.pdf` | Per-gene MAD distribution after median filtering |
| `Final_Filtered_Distribution.pdf` | Expression value histogram after both filters |
| `Soft_Threshold_Diagnostic_Plots.pdf` | Scale-free topology fit and mean connectivity vs. soft-threshold power |
| `Gene_Dendrogram_Module_Colors.pdf` | Gene dendrogram with module colour assignments |

---

## Example

To run this module, place a config file and a Snakefile in the current directory.

### Config (`my_config.yaml`)

**Raw mode** — provide a raw count matrix and let the module normalize it:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    wgcna:
        inputs:
            expression_matrix: "data/counts_matrix.tsv"
        gene_expression:
            input_mode: "raw"
            samples_metadata_path: "data/samples_metadata.tsv"
            failed_qc_path: "data/failed_qc_samples.tsv"   # optional; set to "" to skip
            pathologies:
                - DLBCL
                - FL
```

**Normalized mode** — supply an already-normalized matrix and skip straight to filtering:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    wgcna:
        inputs:
            expression_matrix: "data/normalized_expression.rds"
        gene_expression:
            input_mode: "normalized"
```

### Snakefile

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/wgcna/1.0/config/default.yaml"
configfile: "my_config.yaml"

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/wgcna/1.0/wgcna.smk"

rule all:
    input:
        rules._wgcna_all.input
```

---

## Pipeline overview


```
[raw counts TSV]                  [normalized matrix RDS/TSV]
       │                                      │
       ▼  (raw mode)                          │ (normalized mode)
_wgcna_normalize_expression                   │
  DESeq2 VST + limma batch correction         │
       │                                      │
       ▼                                      │
normalized_expression.tsv                     │
       │                                      │
       └──────────────────┬───────────────────┘
                          ▼
          _wgcna_filter_variance_genes
               median & MAD thresholds
                          │
                          ▼
              filtered_expression.tsv
                          │
                          ▼
          _wgcna_get_coexpression_modules
               soft-threshold selection
               blockwiseModules (WGCNA)
                          │
          ┌───────────────┴────────────────┐
          ▼                                ▼
 coexpression_modules.rds             network.rds
 + coexpression_modules.tsv           (full network object)
 (gene → module colour)
```

---

## Changelog

See the full changelog [here](CHANGELOG.md).