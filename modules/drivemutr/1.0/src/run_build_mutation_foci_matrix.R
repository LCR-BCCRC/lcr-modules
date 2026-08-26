log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tidyr)
library(rlang)

# lambda_suffix() is defined in annotate_lambda.R and used inside
# build_mutation_foci_matrix, so both files must be sourced.
source(file.path(snakemake@scriptdir, "R", "annotate_lambda.R"))
source(file.path(snakemake@scriptdir, "R", "build_mutation_foci_matrix.R"))

gene_data     <- readRDS(snakemake@input[[1]])
lambda_values <- snakemake@params[["lam"]]

# build_mutation_foci_matrix returns a named list of wide data frames
# (one per group type), so annotate_lambda is called with add_df = TRUE.
# The result is a tibble(lambda_name, data) stored as matched_mut_foci_data.
gene_data <- gene_data %>%
  dplyr::mutate(
    matched_mut_foci_data = purrr::map(
      ssm_matched,
      ~ annotate_lambda(df = .x, def = build_mutation_foci_matrix, lambda = lambda_values, add_df = TRUE) %>%
        dplyr::rename(genomic_matrix = data)
    )
  )

saveRDS(gene_data, snakemake@output[[1]])