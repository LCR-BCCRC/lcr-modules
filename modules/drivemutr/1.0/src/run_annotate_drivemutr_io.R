log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(readr)
library(tidyr)

source(file.path(snakemake@scriptdir, "R", "fix_ids.R"))
source(file.path(snakemake@scriptdir, "R", "annotate_drivemutr_io.R"))

ssm <- readr::read_tsv(snakemake@input[["ssm"]], show_col_types = FALSE)
cnv <- readr::read_tsv(snakemake@input[["cnv"]], show_col_types = FALSE)

map_path <- snakemake@params[["sample_id_map_path"]]
sample_id_map <- if (!is.null(map_path) && nzchar(map_path)) {
  readr::read_tsv(map_path, show_col_types = FALSE)
} else {
  NULL
}

target_regions_data <- annotate_drivemutr_io(
  pathology          = snakemake@params[["pathology"]],
  genome_build       = snakemake@params[["genome_build"]],
  genes_regions_list = snakemake@params[["genes_regions_list"]],
  ssm                = ssm,
  cnv                = cnv,
  ashm_regions       = snakemake@params[["ashm_regions"]],
  sample_id_map      = sample_id_map
)

saveRDS(target_regions_data, snakemake@output[[1]])