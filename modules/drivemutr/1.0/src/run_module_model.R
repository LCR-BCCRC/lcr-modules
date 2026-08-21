log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tibble)
library(broom)

source(file.path(snakemake@scriptdir, "R", "module_model.R"))

gene_data <- readRDS(snakemake@input[[1]])

gene_data <- annotate_module_model_results(gene_data)

saveRDS(gene_data, snakemake@output[[1]])