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

gene_data          <- readRDS(snakemake@input[["gene_data"]])
tf_names_expressed <- readr::read_lines(snakemake@input[["tf_names_expressed"]])

gene_data <- annotate_tf_results(gene_data, tf_names_expressed)

saveRDS(gene_data, snakemake@output[[1]])