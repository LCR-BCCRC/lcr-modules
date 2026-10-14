log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(purrr)

source(file.path(snakemake@scriptdir, "R", "get_annotate_cadd_scores.R"))

target_regions_data <- readRDS(snakemake@input[[1]])

cadd_path <- snakemake@params[["cadd_path"]]
cadd_data <- get_CADD_scores(target_regions_data, cadd_path = if (nzchar(cadd_path)) cadd_path else NULL)

target_regions_data <- target_regions_data %>%
  mutate(
    ssm = map(ssm, function(maf) {
      if (is.null(maf) || nrow(maf) == 0) return(maf)
      annotate_CADD_scores(cadd_data, maf)
    })
  )

saveRDS(target_regions_data, snakemake@output[[1]])