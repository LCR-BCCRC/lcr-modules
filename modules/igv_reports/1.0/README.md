**Module Overview**
- **Purpose**: Generate self-contained per-sample HTML variant review reports using [igv-reports](https://github.com/igvteam/igv-reports). Each report embeds tumour and normal BAM pileups and a Gencode gene annotation track for coding variants in bona fide lymphoma driver genes.

**Prerequisites**

The project Snakefile must declare the `reference_files` subworkflow (genomes/{genome_build}/genome_fasta/genome.fa and genomes/{genome_build}/annotations/gencode_annotation-33.gtf are required at runtime).

**Required config**

| Key | Description |
|-----|-------------|
| `inputs.mafs` | Dict mapping a short label to a MAF path pattern. The label becomes the `{tool}` wildcard. |
| `inputs.tumour_bam` | Path pattern for tumour BAM (supports standard sample wildcards) |
| `inputs.normal_bam` | Path pattern for normal BAM |

BAM index files (`.bai`) are inferred by appending `.bai` to the BAM path. Both `.bai` and `.crai` symlinks are created pointing to the same source index so pysam can discover whichever extension it expects.

**Multiple MAF inputs**

Any number of MAF sources can be processed in a single module include by adding entries to `inputs.mafs`:

```yaml
igv_reports:
    inputs:
        mafs:
            sage-only: "results/maf_diff-1.0/99-outputs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage-only.maf"
            slms-3-only: "results/maf_diff-1.0/99-outputs/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.slms-3-only.maf"
```

Each label produces an independent set of HTML reports without needing to include the module twice.

**Outputs** (under `99-outputs/{seq_type}--{genome_build}/`)

One HTML file per tumour/normal pair per tool label:
`{tumour_id}--{normal_id}--{pair_status}.{tool}.html`

Reports are self-contained and can be opened directly in a browser without network access.

**Filtering**

Before creating the report, variants are filtered to coding variants in bona fide driver genes (same gene list and variant class logic as `maf_diff`). This keeps reports focused on actionable variants and avoids extremely large HTML files. The filtered MAF is written to `02-filtered_maf/` and is preserved for inspection.

**Driver gene list**

Defaults to the bundled LLMPP curated tier-1 gene list (`etc/any_tier1_BL_FL_DLBCL_MCL_MZL_PMBL.tsv`). Configurable via `options.driver_genes`, `options.driver_col`, and `options.driver_col_value`.

**Notes**
- The `--fasta` flag embeds reference sequence from the local `reference_files` genome FASTA, avoiding the chromosome naming mismatch that occurs when using `--genome hg19` with grch37 (non-`chr`-prefixed) MAFs and BAMs.
- The Gencode GTF is passed as a track and embedded in the report; gene annotations do not require network access at view time.
- `options.flanking` (default 1000 bp) controls how much sequence context is shown around each variant.
- `options.info_columns` controls which MAF columns appear in the variant table at the top of each report.
