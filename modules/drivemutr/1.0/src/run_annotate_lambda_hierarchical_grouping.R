log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(ggplot2)

source(file.path(snakemake@scriptdir, "R", "annotate_lambda.R"))
source(file.path(snakemake@scriptdir, "R", "hierachical_grouping.R"))
source(file.path(snakemake@scriptdir, "R", "remove_identical_grouped_cols.R"))

gene_data <- readRDS(snakemake@input[[1]])

lambda_values <- snakemake@params[["lam"]]

gene_data <- gene_data %>%
  mutate(
    results     = purrr::map(ssm, ~ annotate_lambda(df = .x, def = hierachical_grouping, lambda = lambda_values)),
    ssm         = purrr::map(results, "maf"),
    height_plot = purrr::map(results, "plot")
  ) %>%
  dplyr::select(-results)

if (length(lambda_values) > 1) {
  gene_data <- gene_data %>%
    mutate(ssm = purrr::map(ssm, remove_identical_grouped_cols))
}

plots <- purrr::compact(gene_data$height_plot)
pdf(snakemake@output[["height_plot_pdf"]], width = 10, height = 6)
for (p in plots) tryCatch(print(p), error = function(e) NULL)
dev.off()

saveRDS(gene_data, snakemake@output[["rds"]])