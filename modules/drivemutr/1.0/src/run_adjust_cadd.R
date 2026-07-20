log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)

source(file.path(snakemake@scriptdir, "R", "annotate_lambda.R"))
source(file.path(snakemake@scriptdir, "R", "adjust_cadd.R"))

gene_data <- readRDS(snakemake@input[[1]])

lambda_values <- snakemake@params[["lam"]]

gene_data <- gene_data %>%
  dplyr::mutate(
    ssm = purrr::map(ssm, ~ annotate_lambda(df = .x, def = adjust_cadd, lambda = lambda_values))
  )

saveRDS(gene_data, snakemake@output[[1]])