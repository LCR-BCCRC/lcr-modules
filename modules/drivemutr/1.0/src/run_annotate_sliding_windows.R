log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(rlang)

source(file.path(snakemake@scriptdir, "R", "annotate_sliding_windows.R"))

gene_data <- readRDS(snakemake@input[[1]])

gene_data <- gene_data %>%
  mutate(
    ssm = purrr::map(ssm, ~ split_region_custom_col(
      ssm_df              = .x,
      window_size         = snakemake@params[["window_size"]],
      max_range           = snakemake@params[["max_range"]],
      slide_by            = snakemake@params[["slide_by"]],
      min_mutations_count = snakemake@params[["min_mutations_count"]],
      min_group_range     = snakemake@params[["min_group_range"]]
    ))
  )

saveRDS(gene_data, snakemake@output[[1]])