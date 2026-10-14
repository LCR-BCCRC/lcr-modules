log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(readr)

source(file.path(snakemake@scriptdir, "R", "fix_ids.R"))
source(file.path(snakemake@scriptdir, "R", "assign_module_color_and_expression.R"))

gene_data              <- readRDS(snakemake@input[["gene_data"]])
mergedColors           <- readRDS(snakemake@input[["coexpression_modules"]])
norm <- readr::read_tsv(snakemake@input[["normalized_expression"]], show_col_types = FALSE)
if (!"hgnc_symbol" %in% colnames(norm)) {
  stop("normalized_expression is missing the 'hgnc_symbol' column expected from the WGCNA module")
}
sample_cols    <- setdiff(colnames(norm), c("gene_id", "hgnc_symbol"))
gene_by_sample <- as.matrix(norm[, sample_cols, drop = FALSE])
rownames(gene_by_sample) <- norm[["hgnc_symbol"]]
expression_matrix <- t(gene_by_sample)  # samples x genes

map_path <- snakemake@params[["sample_id_map_path"]]
sample_id_aliases <- if (!is.null(map_path) && nzchar(map_path)) {
  readr::read_tsv(map_path, show_col_types = FALSE) %>%
    tidyr::pivot_longer(
      cols      = c(mrna_sample_id, genome_sample_id, capture_sample_id),
      names_to  = "id_type",
      values_to = "all_ids"
    )
} else {
  NULL
}

gene_data <- assign_module_color_and_expression(
  target_regions_df          = gene_data,
  mergedColors               = mergedColors,
  expression_matrix          = expression_matrix,
  sample_id_aliases          = sample_id_aliases
)

saveRDS(gene_data, snakemake@output[[1]])