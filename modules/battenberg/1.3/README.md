**Module Overview**
- **Purpose**: Two-stage Battenberg WGS module: heavy preprocessing produces reusable `.tab` files, then Snakemake runs per-ploidy fitting jobs using those tables.

**Config: `options.ploidy_runs`**
- **Type**: list of strings.
- **Default**: `["1.6-4.8"]` (matches Battenberg defaults).
- **Format**: each entry is `MIN-MAX` (for example `"1.6-4.8"`).
- **Behavior**: Snakemake expands the pipeline to run one lightweight fit job per list entry; the first list entry is treated as the canonical ploidy for downstream outputs.

**Runner CLI: `--ploidy_constraint`**
- **Usage**: pass a constraint like `--ploidy_constraint 1.6-4.8` to the R runner (`src/battenberg_wgs_hg38.R`).
- **Effect**: the runner parses the `MIN-MAX` string and sets `min_ploidy`/`max_ploidy` for that fit.
- **Backward compatibility**: the runner still accepts `--min_ploidy` and `--max_ploidy`; explicit `--min_ploidy`/`--max_ploidy` override a `--ploidy_constraint` value.

**Notes**
- **Preprocessing**: heavy preprocessing outputs (allele counts, impute outputs, `.tab` files) are written into a `preprocess/` subdirectory and intentionally preserved so they can be reused by multiple fit runs.
- **Cleanup**: `_battenberg_cleanup` removes the preserved `preprocess/*.tab` files only after all per-ploidy fit jobs complete.
- **Canonical outputs**: downstream rules (plots, seg files, etc.) reference the canonical ploidy defined as the first entry in `options.ploidy_runs`.

**Examples**
- **Snakemake (dry-run)**: `snakemake -n -s modules/battenberg/1.2/battenberg.smk _battenberg_all`
- **Direct R runner (single fit, skipping preprocessing)**:
  `Rscript modules/battenberg/1.2/src/battenberg_wgs_hg38.R -t TUMOUR -n NORMAL --tb tumour.bam --nb normal.bam -o outdir --ploidy_constraint 1.6-4.8 --skip_preprocessing --skip_allelecount`

**Multi-run example**
- To run two fits per input (e.g. the default plus a broader range), set in `config/default.yaml`:

```yaml
options:
  ploidy_runs: ["1.6-4.8", "2.0-5.0"]
```

- Snakemake will expand the DAG to run one preprocess step, then two lightweight fits per pair producing outputs under `ploidy_1.6-4.8/` and `ploidy_2.0-5.0/`. The first entry (`1.6-4.8`) remains the canonical run used for downstream default outputs.

If you want a different canonical ploidy than the first list item, change the order of `options.ploidy_runs` in the module config.

**Canonical ploidy selection**
- The module treats the first entry of `options.ploidy_runs` as the canonical ploidy. Downstream outputs (plots, default seg/txt exports) use the canonical run. To change which fit is canonical, reorder the `ploidy_runs` list in `config/default.yaml`.
