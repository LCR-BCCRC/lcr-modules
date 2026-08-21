log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

source(file.path(snakemake@scriptdir, "R", "annotate_sanity_check_plot.R"))

gene_data <- readRDS(snakemake@input[[1]])

gene_data <- annotate_sanity_check_plot(gene_data)

all_plots <- purrr::compact(purrr::flatten(
  purrr::map(gene_data$matched_mut_foci_data, function(mfd) {
    if (is.null(mfd) || !"Sanity_Check_Plot" %in% names(mfd)) return(list())
    as.list(mfd$Sanity_Check_Plot)
  })
))
pdf(snakemake@output[["sanity_check_pdf"]], width = 8, height = 6)
for (p in all_plots) tryCatch(print(p), error = function(e) NULL)
dev.off()

saveRDS(gene_data, snakemake@output[["rds"]])