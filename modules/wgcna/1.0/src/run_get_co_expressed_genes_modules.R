log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(con) })

library(dplyr)
library(tibble)
library(readr)
library(WGCNA)

source(file.path(snakemake@params[["script_dir"]], "R", "get_co_expressed_genes_modules.R"))

dir.create(dirname(snakemake@output[["modules_rds"]]), recursive = TRUE, showWarnings = FALSE)

WGCNA::allowWGCNAThreads(nThreads = snakemake@threads)

expression_matrix <- read_tsv(snakemake@input[["expression"]], show_col_types = FALSE) %>%
  column_to_rownames("sample_id") %>%
  as.matrix()

result <- get_co_expressed_genes_modules(
  expression_matrix = expression_matrix,
  powers            = as.integer(snakemake@params[["powers"]]),
  cor_method        = snakemake@params[["cor_method"]],
  network_type      = snakemake@params[["network_type"]],
  plots_dir         = snakemake@params[["plots_dir"]]
)

saveRDS(result$module_colors, snakemake@output[["modules_rds"]])
saveRDS(result$network,       snakemake@output[["network_rds"]])

# TSV: one row per gene with its assigned module colour
modules_tsv <- data.frame(
  gene   = names(result$module_colors),
  module = result$module_colors,
  row.names = NULL
)
write_tsv(modules_tsv, snakemake@output[["modules_tsv"]])