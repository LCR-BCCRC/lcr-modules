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
high_var_expression_matrix <- readr::read_tsv(snakemake@input[["filtered_expression"]], show_col_types = FALSE) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()

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
  high_var_expression_matrix = high_var_expression_matrix,
  sample_id_aliases          = sample_id_aliases
)

saveRDS(gene_data, snakemake@output[[1]])