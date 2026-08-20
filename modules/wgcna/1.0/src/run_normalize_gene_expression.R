log_file <- snakemake@log[[1]]
con <- file(log_file, open = "wt")
sink(con, type = "output")
sink(con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(con) })

library(readr)
library(dplyr)
library(tibble)
library(DESeq2)
library(limma)
library(BiocParallel)
library(GAMBLR.helpers)

source(file.path(snakemake@params[["script_dir"]], "R", "normalize_gene_expression.R"))

dir.create(dirname(snakemake@output[["tsv"]]), recursive = TRUE, showWarnings = FALSE)

mat <- normalize_gene_expression_by_pathology(
  raw_gene_exp_path     = snakemake@params[["raw_matrix_path"]],
  failed_qc_path        = snakemake@params[["failed_qc_path"]],
  samples_metadata_path = snakemake@params[["samples_metadata_path"]],
  pathology_arg         = snakemake@params[["pathology"]],
  plots_dir             = snakemake@params[["plots_dir"]],
  batches               = snakemake@params[["batches"]],
  cohort_var            = snakemake@params[["cohort_var"]],
  bio_var               = snakemake@params[["bio_var"]],
  min_reads             = snakemake@params[["min_reads"]],
  threads               = snakemake@threads
)

# mat is a genes × samples data frame with hgnc_symbol column
write_tsv(mat, snakemake@output[["tsv"]])
