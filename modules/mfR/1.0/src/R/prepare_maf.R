#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# prepare_maf.R  (run under Snakemake `script:` directive)
# Subset a (multi-seq_type) master MAF to one sample_set's samples and keep
# only non-coding mutations, by dropping the configured coding
# Variant_Classification values. Mirrors anti_join(all_ssm, coding_ssm).
# ---------------------------------------------------------------------------

log <- file(snakemake@log[["log"]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

maf_paths        <- snakemake@input[["maf"]]            # one per seq_type
sample_sets_path <- snakemake@input[["sample_sets"]]
sample_set       <- snakemake@wildcards[["sample_set"]]
out_maf          <- snakemake@output[["maf"]]

coding_classes   <- snakemake@params[["coding_classes"]]     # character vector (maybe length 0)
sid_col          <- snakemake@params[["sample_id_col"]]
ss_col           <- snakemake@params[["sample_set_col"]]
maf_sample_col   <- snakemake@params[["maf_sample_col"]]

# Samples belonging to this sample_set
ss <- read_tsv(sample_sets_path, show_col_types = FALSE)
stopifnot(all(c(sid_col, ss_col) %in% names(ss)))
ids <- ss[[sid_col]][ ss[[ss_col]] == sample_set ]
message(sprintf("sample_set '%s': %d samples", sample_set, length(ids)))

# Combine seq_type MAFs, restrict to this sample_set
maf <- map_dfr(as.character(maf_paths), ~ read_tsv(.x, show_col_types = FALSE))
stopifnot(maf_sample_col %in% names(maf))
maf <- maf %>% filter(.data[[maf_sample_col]] %in% ids)

# Keep non-coding: drop the configured coding Variant_Classification values
if (length(coding_classes) > 0 && "Variant_Classification" %in% names(maf)) {
  before <- nrow(maf)
  maf <- maf %>% filter(!Variant_Classification %in% coding_classes)
  message(sprintf("Dropped %d coding rows; %d non-coding remain",
                  before - nrow(maf), nrow(maf)))
}

dir.create(dirname(out_maf), recursive = TRUE, showWarnings = FALSE)
write_tsv(maf, out_maf)
message(sprintf("Wrote %s (%d rows)", out_maf, nrow(maf)))
