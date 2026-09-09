log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tibble)

source(file.path(snakemake@scriptdir, "R", "merge_multi_ashm_regions_gene.R"))

annotated_data <- readRDS(snakemake@input[[1]])

merged_data <- merge_multi_ashm_regions_gene(annotated_data)

out_dir <- snakemake@output[[1]]

for (i in seq_len(nrow(merged_data))) {
  gene     <- merged_data$Hugo_Symbol[i]
  gene_dir <- file.path(out_dir, gene)
  dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(merged_data[i, ], file.path(gene_dir, "data.rds"))
}

message("Split ", nrow(merged_data), " genes into ", out_dir)