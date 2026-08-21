log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifbreakR)
library(MotifDb)
library(BiocParallel)

source(file.path(snakemake@scriptdir, "R", "get_tf_results.R"))

gene_data            <- readRDS(snakemake@input[["gene_data"]])
filtered_expression  <- readr::read_tsv(snakemake@input[["filtered_expression"]], show_col_types = FALSE) %>%
  tibble::column_to_rownames("sample_id") %>%
  as.matrix()
tf_names             <- readr::read_lines(snakemake@input[["tf_names_file"]])
expression_threshold <- snakemake@params[["expression_threshold"]]

# Derive expressed TF names: keep TFs present in the expression matrix whose
# mean VST value meets the threshold
tf_exp    <- filtered_expression[, colnames(filtered_expression) %in% tf_names, drop = FALSE]
tf_means  <- colMeans(tf_exp, na.rm = TRUE)
tf_names_expressed <- names(tf_means)[tf_means >= expression_threshold]

gene_data <- annotate_tf_results(gene_data, tf_names_expressed)

saveRDS(gene_data, snakemake@output[[1]])