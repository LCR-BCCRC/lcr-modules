log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(pre)
library(iml)

source(file.path(snakemake@scriptdir, "R", "rulefit_shap_model.R"))

gene_data <- readRDS(snakemake@input[[1]])

gene_data <- get_rulefit_shapley_results(gene_data)

all_plots <- purrr::compact(purrr::flatten(
  purrr::map(gene_data$matched_mut_foci_data, function(mfd) {
    if (is.null(mfd) || !"Interaction_Results" %in% names(mfd)) return(list())
    purrr::flatten(
      purrr::map(mfd$Interaction_Results, function(ir) {
        if (is.null(ir) || nrow(ir) == 0) return(list())
        c(ir$SHAP_Mean_Plot, ir$SHAP_Beeswarm_Plot)
      })
    )
  })
))
pdf(snakemake@output[["shap_pdf"]], width = 8, height = 6)
for (p in all_plots) tryCatch(print(p), error = function(e) NULL)
dev.off()

saveRDS(gene_data, snakemake@output[["rds"]])