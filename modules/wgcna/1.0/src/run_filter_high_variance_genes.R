log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(con) })

library(ggplot2)
library(readr)
library(dplyr)
library(tibble)

source(file.path(snakemake@params[["script_dir"]], "lib", "filter_high_variance_genes.R"))

dir.create(dirname(snakemake@output[["tsv"]]), recursive = TRUE, showWarnings = FALSE)

input_file <- snakemake@input[["expression"]]
if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
  expression_matrix <- readRDS(input_file)
} else {
  expression_matrix <- read_tsv(input_file, show_col_types = FALSE)
}

filtered <- filter_high_variance_genes(
  expression_matrix = expression_matrix,
  median_threshold  = snakemake@params[["median_threshold"]],
  mad_threshold     = snakemake@params[["mad_threshold"]],
  plots_dir         = snakemake@params[["plots_dir"]]
)

filtered_tsv <- tibble::rownames_to_column(as.data.frame(filtered), var = "sample_id")
write_tsv(filtered_tsv, snakemake@output[["tsv"]])