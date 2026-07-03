log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tidyr)
library(readr)

source(file.path(snakemake@scriptdir, "R", "get_ucsc_custom_tracks.R"))

# Load every per-gene sanity_check_plot.rds and bind into one tibble.
# snakemake@input is a list of paths (one per gene, resolved by per_gene_data).
target_results <- purrr::map(unlist(snakemake@input), readRDS) %>%
  dplyr::bind_rows()

outdir <- snakemake@params[["outdir"]]

get_ucsc_custom_tracks(target_results, outdir)