log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tidyr)

source(file.path(snakemake@scriptdir, "R", "filter_matched_DNA_RNA_samples.R"))

gene_data <- readRDS(snakemake@input[[1]])

gene_data <- filter_matched_DNA_RNA_samples(gene_data)

saveRDS(gene_data, snakemake@output[[1]])