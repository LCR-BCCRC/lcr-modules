#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cluster_foci.R  (run under Snakemake `script:` directive)
# Hierarchical clustering of non-coding mutation positions into "foci" for a
# SINGLE chromosome of a SINGLE sample_set. Port of the notebook's
# `hierachical_grouping_silhouette()`.
#
# The input MAF is already scoped to this sample_set x chromosome by the
# upstream _mfr_extract_chrom rule (tabix-based extraction), so unlike the
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
#   * dist() is O(n^2) in the number of unique positions, which is the main OOM
#     risk once a sample_set pools ~700 WGS samples onto one chromosome. Before
#     clustering, positions are split into gap-delimited chunks wherever two
#     consecutive unique positions are more than h_max apart, and each chunk is
#     clustered independently -- see split_into_gap_chunks() below for why that
#     is exact (not approximate) for cut heights <= h_max. This bounds dist()
#     to the largest chunk rather than the whole chromosome, but is not a hard
#     memory ceiling: a chromosome region where consecutive positions across
#     all samples never exceed h_max apart (e.g. a shared hotspot) still forms
#     one large chunk. Watch the per-chunk sizes this script logs; if the
#     largest chunk is still too big, h_max itself needs to come down.
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

# ---- clustering for one gap-chunk of positions -----------------------------
# Identical math to the pre-chunking version of this script: builds one
# dist()/hclust() over this chunk's unique positions and picks the cut height
# maximizing mean silhouette width. Returns chunk-LOCAL group ids (starting at
# 1) -- callers must offset them to stay unique across chunks.
cluster_one_chunk <- function(maf, pos_col, dist_method, hclust_method,
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
    # silhouette() needs >1 cluster AND at least one non-singleton cluster --
    # with every point in its own cluster (common at h_min before any merge
    # has happened) it returns a bare NA instead of a matrix, and indexing
    # that with [, "sil_width"] throws "incorrect number of dimensions".
    if (length(unique(lab)) > 1 && length(unique(lab)) < length(lab)) {
      # Belt-and-suspenders: the condition above rules out the two known
      # degenerate cases, but wrap anyway so any other silhouette() edge
      # case can't abort the job either (matches the cutree() tryCatch above).
      sil_means[i] <- tryCatch(
        mean(silhouette(lab, d)[, "sil_width"]),
        error = function(e) NA_real_
      )
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

# ---- split unique positions into gap-delimited chunks ----------------------
# A chunk boundary is placed between two consecutive unique positions whenever
# they are more than h_max apart. This is exact, not approximate: for any two
# point sets that lie entirely on either side of such a gap G, every pairwise
# distance between a point on the left and a point on the right is >= G (by
# construction of 1D sorted coordinates), so the merge height for single,
# complete, average or centroid linkage -- and the centroid-distance term
# ward.D2 is built from -- between any left-subset and any right-subset is
# also >= G. Since G > h_max, no cut at h <= h_max can ever join points from
# different chunks, so clustering each chunk independently is equivalent to
# clustering the whole chromosome at any h <= h_max. Returns a chunk id per
# element of `positions` (already assumed sorted, unique).
split_into_gap_chunks <- function(positions, h_max) {
  n <- length(positions)
  if (n < 2) return(rep(1L, n))
  gaps <- diff(positions)
  cumsum(c(TRUE, gaps > h_max))
}

# ---- cluster a whole chromosome by clustering each gap-chunk separately ----
cluster_chromosome <- function(maf, pos_col, dist_method, hclust_method,
                               h_min, h_max) {
  if (nrow(maf) == 0) {
    maf$group <- integer(0)
    return(list(maf = maf,
                plot_df = data.frame(h = integer(0), avg = numeric(0),
                                     chunk = integer(0), best_h = integer(0)),
                n_chunks = 0L))
  }

  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  chunk_of_position <- split_into_gap_chunks(positions, h_max)
  # Numeric index match (not a character-keyed lookup): positions here are
  # exactly the values maf[[pos_col]] was derived from, so no risk of a
  # scientific-notation mismatch on large genomic coordinates.
  maf$.chunk <- chunk_of_position[match(maf[[pos_col]], positions)]

  n_chunks <- max(chunk_of_position)
  message(sprintf(
    "  %d unique position(s) split into %d gap-chunk(s) (h_max = %d); largest chunk = %d position(s)",
    length(positions), n_chunks, h_max, max(table(chunk_of_position))
  ))

  maf_chunks   <- vector("list", n_chunks)
  plot_chunks  <- vector("list", n_chunks)
  group_offset <- 0L

  for (cid in seq_len(n_chunks)) {
    sub <- maf[maf$.chunk == cid, setdiff(names(maf), ".chunk"), drop = FALSE]
    res <- cluster_one_chunk(sub, pos_col, dist_method, hclust_method, h_min, h_max)

    local_max <- if (nrow(res$maf) > 0) max(res$maf$group) else 0L
    res$maf$group <- res$maf$group + group_offset
    maf_chunks[[cid]] <- res$maf

    if (nrow(res$plot_df) > 0) {
      plot_chunks[[cid]] <- cbind(res$plot_df, chunk = cid, best_h = res$best_h)
    }

    group_offset <- group_offset + local_max
  }

  plot_df <- bind_rows(plot_chunks)
  if (nrow(plot_df) == 0) {
    plot_df <- data.frame(h = integer(0), avg = numeric(0),
                          chunk = integer(0), best_h = integer(0))
  }

  list(maf = bind_rows(maf_chunks), plot_df = plot_df, n_chunks = n_chunks)
}

# ---- run -------------------------------------------------------------------
maf_chr <- read_tsv(maf_path, show_col_types = FALSE, progress = FALSE)
message(sprintf("chr %s: %d rows, %d unique positions",
                chrom, nrow(maf_chr), length(unique(maf_chr[[pos_col]]))))

res <- cluster_chromosome(maf_chr, pos_col, dist_method,
                          hclust_method, h_min, h_max)

write_tsv(res$maf, out_tsv, na = "")

p <- ggplot(res$plot_df, aes(x = h, y = avg)) +
  geom_line() +
  geom_point(shape = 16) +
  labs(x = "Height (h)", y = "Mean silhouette width",
       title = sprintf("Foci clustering quality — chr %s (%d gap-chunk(s))",
                       chrom, res$n_chunks)) +
  theme_minimal()
if (res$n_chunks > 1) {
  p <- p + facet_wrap(~ chunk, scales = "free_y")
}
ggsave(out_plot, p, width = 7, height = 4)

message(sprintf("chr %s: %d gap-chunk(s), %d foci total -> %s",
                chrom, res$n_chunks, length(unique(res$maf$group)), out_tsv))
