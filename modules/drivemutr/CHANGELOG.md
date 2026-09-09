# Changelog

All notable changes to the `drivemutr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0] - 2026-07-03

This release was authored by Houman Layegh Mirhosseini.

- Initial release of the `drivemutr` module: a Level 4, cohort-level, per-gene workflow that searches aSHM (aberrant Somatic HyperMutation) and user-supplied target regions for non-coding regulatory driver mutation foci — clusters of mutations whose presence is associated with a change in the expression of the host gene.
- The module consumes a single cohort-level MAF (`ssm_maf`) and copy-number matrix (`cnv_matrix`), together with the co-expression module assignments (`wgcna_coexpression_modules`) and filtered expression matrix (`wgcna_filtered_expression`) from the `WGCNA` module, plus a TF gene-symbol list (`tf_names_file`). Optional inputs are a DNA↔RNA linkage table (`sample_id_map`) and a directory of local CADD bigwig files (`cadd_dir`); when `cadd_dir` is null, CADD scores are fetched live from UCSC. `genome_build`, `pathology`, and the target-region definition (`ashm_regions` / `genes_regions_list`) are supplied through `reference_params`.
- Because the module is cohort-level and per-gene rather than per-sample, it does not use tumour/normal pairing and no `pairing_config` is required. Its working wildcards are `{gene}` and `{lam}` (one hierarchical-grouping lambda value) rather than the usual `seq_type`/`sample_id`. Target genes are scattered via a Snakemake `checkpoint` (`_drivemutr_merge_and_split_genes`) and gathered back into the aggregated tracks.
- Per gene, the pipeline scores mutations for deleteriousness with CADD; groups nearby positions into foci by hierarchical clustering (controlled by `lambda`), optionally refined with a sliding window (`custom_sliding_window`); links each mutated sample to its WGCNA co-expression module and host-gene expression, keeping only samples with matched DNA and RNA; fits module-eigengene, genomic, and RuleFit + Shapley interaction models to find foci significantly associated with expression; and runs motifbreakR TF-motif disruption analysis with a Fisher exact test for enrichment of expressed transcription factors.
- Outputs are UCSC custom track TSVs (`Mutation_Points`, `Mutation_Blocks`, `Transcription_Factors`), one set per `lambda` value, symlinked under `99-outputs/ucsc_custom_tracks/`. Per-gene diagnostic plots (`height_plot.pdf`, `shap_plots.pdf`, `sanity_check_plot.pdf`) and the terminal `final_results.rds` are retained under the per-gene working directory; the intermediate per-gene `.rds` files (one per pipeline step) are marked `temp()` and removed automatically once the gene's `final_results.rds` is produced.
- All rules run in a dedicated conda environment (`envs/drivemutr-1.0.yaml`) or, with `--software-deployment-method apptainer`, the module container (`docker://ghcr.io/LCR-BCCRC/lcr-scripts/drivemutr:1.0`).