#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cluster_foci.R  (run under Snakemake `script:` directive)
# Hierarchical clustering of non-coding mutation positions into "foci" for a
# SINGLE chromosome of a SINGLE sample_set. Port of the notebook's
# `hierachical_grouping_silhouette()`.
#
# The input MAF is already scoped to this sample_set x chromosome by the
# upstream _mfR_extract_chrom rule (tabix-based extraction), so unlike the
# old mutation_foci module this script no longer reads a genome-wide or
# whole-sample_set file and filters it down -- it just reads the (small)
# already-scoped file directly.
#
# Outputs: (1) the MAF rows for this chromosome with an added `group` (focus)
#              column; (2) a silhouette-vs-height diagnostic plot.
#
# Notes:
#   * Clusters are built on UNIQUE positions; every mutation at a position
#     inherits that position's focus id (foci are positional).
#   * centroid/median linkage can yield non-monotonic merge heights, which makes
#     height-based cutting unreliable; each cut is wrapped in tryCatch so one bad
#     height (or chromosome) cannot abort the job. Switch hclust_method to
#     "average"/"ward.D2" in the config if foci look wrong.
#   * dist() is O(n^2) in the number of unique positions; see CHANGELOG for the
#     gap-based pre-blocking optimisation if a dense chromosome OOMs.
# ---------------------------------------------------------------------------

log <- file(snakemake@log[["log"]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(cluster)    # silhouette()
  library(ggplot2)
})

maf_path      <- snakemake@input[["maf"]]
chrom         <- snakemake@wildcards[["chrom"]]
out_tsv       <- snakemake@output[["tsv"]]
out_plot      <- snakemake@output[["plot"]]

pos_col       <- snakemake@params[["pos_col"]]
dist_method   <- snakemake@params[["dist_method"]]
hclust_method <- snakemake@params[["hclust_method"]]
h_min         <- as.integer(snakemake@params[["h_min"]])
h_max         <- as.integer(snakemake@params[["h_max"]])

# ---- clustering for one chromosome ----------------------------------------
cluster_one_chromosome <- function(maf, pos_col, dist_method, hclust_method,
                                   h_min, h_max) {
  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  n <- length(positions)

  # Degenerate: 0 or 1 unique position -> single focus, no silhouette.
  if (n < 2) {
    maf$group <- if (n == 0) integer(0) else 1L
    return(list(maf = maf,
                plot_df = data.frame(h = integer(0), avg = numeric(0)),
                best_h = NA_integer_))
  }

  d  <- dist(positions, method = dist_method)
  hc <- hclust(d, method = hclust_method)

  h_values  <- seq.int(h_min, h_max)
  sil_means <- rep(NA_real_, length(h_values))

  for (i in seq_along(h_values)) {
    lab <- tryCatch(cutree(hc, h = h_values[i]), error = function(e) NULL)
    if (is.null(lab)) next
    if (length(unique(lab)) > 1) {                 # silhouette needs > 1 cluster
      sil_means[i] <- mean(silhouette(lab, d)[, "sil_width"])
    }
  }

  best_h <- if (all(is.na(sil_means))) h_min else h_values[which.max(sil_means)]

  labels <- tryCatch(cutree(hc, h = best_h),
                     error = function(e) seq_along(positions))   # fallback: singletons
  lab_df <- setNames(data.frame(positions, labels), c(pos_col, "group"))
  maf <- left_join(maf, lab_df, by = pos_col)

  list(maf = maf,
       plot_df = data.frame(h = h_values, avg = sil_means),
       best_h = best_h)
}

# ---- run -------------------------------------------------------------------
maf_chr <- read_tsv(maf_path, show_col_types = FALSE, progress = FALSE)
message(sprintf("chr %s: %d rows, %d unique positions",
                chrom, nrow(maf_chr), length(unique(maf_chr[[pos_col]]))))

res <- cluster_one_chromosome(maf_chr, pos_col, dist_method,
                              hclust_method, h_min, h_max)

write_tsv(res$maf, out_tsv, na = "")

p <- ggplot(res$plot_df, aes(x = h, y = avg)) +
  geom_line() +
  geom_point(shape = 16) +
  labs(x = "Height (h)", y = "Mean silhouette width",
       title = sprintf("Foci clustering quality — chr %s", chrom)) +
  theme_minimal()
ggsave(out_plot, p, width = 7, height = 4)

message(sprintf("best_h = %s -> %s", as.character(res$best_h), out_tsv))
