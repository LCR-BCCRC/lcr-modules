#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# cluster_foci.R  (run under Snakemake `script:` directive)
# Hierarchical clustering of non-coding mutation positions into "foci" for a
# SINGLE chromosome of a SINGLE sample_set. Supports three interchangeable
# clustering methods, selected via the `clustering_method` config option:
#   * "silhouette" -- pick one cut height per gap-chunk maximizing mean
#     silhouette width (port of the notebook's hierachical_grouping_silhouette()).
#   * "density"    -- pick one cut height per gap-chunk maximizing a
#     density/fragmentation-penalized score (port of hierachical_grouping()'s
#     avg_dist score; math unchanged from Houman's original).
#   * "dynamic"    -- dynamicTreeCut::cutreeDynamic() resolves cluster
#     boundaries per branch, so different foci in the same chunk can end up
#     with different effective heights, instead of one height for the chunk.
#
# The input MAF is already scoped to this sample_set x chromosome by the
# upstream _mfr_extract_chrom rule (tabix-based extraction).
#
# Outputs: (1) the MAF rows for this chromosome with an added `group` (focus)
#              column; (2) h_distribution.pdf, a per-chromosome diagnostic
#              plot -- for silhouette/density, a histogram of the chosen
#              best_h per gap-chunk; for dynamic, a histogram of the
#              effective cut height per FOCUS (no single best_h per chunk
#              exists for that method). (3) h_score.pdf -- for
#              silhouette/density, one point per chunk (its best_h against
#              the score achieved there); for dynamic, a placeholder (this
#              method has no native score to plot -- see notes below).
#              Both plots are skipped (cheap placeholders written instead,
#              so the declared outputs still exist) when make_plots is False.
#
# gap vs. h_max / cut_height -- IMPORTANT:
#   `gap` (chunk-splitting boundary) is a SEPARATE parameter from the height
#   ceiling each method actually searches/considers (h_max for
#   silhouette/density, cut_height for dynamic). For silhouette/density,
#   gap-based chunking is EXACT (not approximate) for any gap >= h_max --
#   see split_into_gap_chunks() below for the proof. That means there is NO
#   benefit to setting gap higher than h_max: results are identical either
#   way, but a bigger gap makes chunks bigger, and dist()/hclust() cost is
#   quadratic in chunk size, so gap should be set close to h_max, not
#   arbitrarily above it. (Caveat: with hclust_method = "ward.D2", the exact
#   bound is gap >= h_max * sqrt(2), not just >= h_max, due to how that
#   linkage's merge-height scaling behaves for small clusters.)
#   For the dynamic method, this exactness proof does NOT transfer --
#   cutreeDynamic has no fixed height ceiling of its own the way a
#   single-height cutree() does, so gap-based chunking there is a pragmatic
#   memory bound, not a proven-equivalent optimization. `cut_height` is
#   still passed to cutreeDynamic() as its own `cutHeight` argument, so each
#   chunk's own search stays bounded independent of what `gap` is set to.
#   Whichever method is active, a warning is printed if gap is set below the
#   method's own height ceiling, since that silently breaks the guarantees
#   above.
#
# Notes:
#   * Clusters are built on UNIQUE positions; every mutation at a position
#     inherits that position's focus id (foci are positional).
#   * centroid/median linkage can yield non-monotonic merge heights. For
#     silhouette/density this only matters for the single tryCatch-guarded
#     cutree() call at best_h, which is safe either way. For the dynamic
#     method, cutreeDynamic's branch-jump heuristics assume roughly
#     monotonic merge heights, so "average" or "ward.D2" is a safer choice
#     of hclust_method than "centroid" when clustering_method == "dynamic".
#   * dist()/hclust() is O(n^2) in the number of unique positions per chunk,
#     the main OOM risk once a sample_set pools many WGS samples onto one
#     chromosome; the density method additionally materializes as.matrix(d)
#     (roughly 2x the memory of the dist object alone) and repeats an O(k^2)
#     submatrix extraction per cluster on every tested height. Watch the
#     per-chunk sizes this script logs; if the largest chunk is too big,
#     `gap` needs to come down -- NOT up, for the reasons above.
#   * cutree()-based methods (silhouette/density) renumber each chunk's
#     cluster ids 1..k by ascending minimum position after clustering, since
#     cutree()'s raw label order is not itself guaranteed to follow spatial
#     order. cutreeDynamic's raw label order reflects the order it
#     discovers/finalizes branches while walking the dendrogram, which is
#     even less likely to follow spatial order (especially with the hybrid
#     method's PAM reassignment step, which can pull individual points into
#     a non-adjacent cluster) -- the same renumbering is applied there too.
#   * Dynamic method only, optional: max_abs_core_scatter/min_abs_gap bound
#     cluster tightness in actual bp rather than by point count the way
#     min_cluster_size does -- e.g. two mutations 800bp apart in an
#     otherwise-empty chunk currently form a valid 2-point focus purely
#     because min_cluster_size = 2 allows any pair, regardless of distance.
#     Setting max_abs_core_scatter rejects that kind of isolated, unbridged
#     pair while still allowing a real chain of closely-spaced mutations to
#     span a similar total distance. Left unset (NULL) by default, in which
#     case cutreeDynamic computes its own internal default for either one.
# ---------------------------------------------------------------------------

log <- file(snakemake@log[["log"]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(cluster)          # silhouette()
  library(dynamicTreeCut)   # cutreeDynamic()
  library(ggplot2)
})

maf_path       <- snakemake@input[["maf"]]
chrom          <- snakemake@wildcards[["chrom"]]
out_tsv        <- snakemake@output[["tsv"]]
out_plot       <- snakemake@output[["plot"]]
out_score_plot <- snakemake@output[["score_plot"]]

pos_col           <- snakemake@params[["pos_col"]]
dist_method       <- snakemake@params[["dist_method"]]
hclust_method     <- snakemake@params[["hclust_method"]]
clustering_method <- snakemake@params[["clustering_method"]]
gap               <- as.integer(snakemake@params[["gap"]])
h_min             <- as.integer(snakemake@params[["h_min"]])
h_max             <- as.integer(snakemake@params[["h_max"]])
cut_height        <- as.integer(snakemake@params[["cut_height"]])
min_cluster_size  <- as.integer(snakemake@params[["min_cluster_size"]])
deep_split        <- as.integer(snakemake@params[["deep_split"]])
make_plots        <- as.logical(snakemake@params[["make_plots"]])

# Optional (NULL/unset by default -- cutreeDynamic computes its own sensible
# default for either one internally when omitted). Read defensively: as.numeric(NULL)
# silently returns a length-0 vector rather than NA, so is.null() must be
# checked first or downstream `if (!is.null(...))` guards would misbehave.
.read_optional_numeric <- function(x) if (is.null(x)) NULL else as.numeric(x)
max_abs_core_scatter <- .read_optional_numeric(snakemake@params[["max_abs_core_scatter"]])
min_abs_gap          <- .read_optional_numeric(snakemake@params[["min_abs_gap"]])

if (!clustering_method %in% c("silhouette", "density", "dynamic")) {
  stop(sprintf(
    "clustering_method must be one of \"silhouette\", \"density\", \"dynamic\" (got \"%s\")",
    clustering_method
  ))
}

# See header note: gap must be >= the active method's own height ceiling to
# keep gap-based chunking from silently changing results (silhouette/
# density) or over-constraining each chunk's search (dynamic).
.height_ceiling <- if (clustering_method == "dynamic") cut_height else h_max
if (gap < .height_ceiling) {
  warning(sprintf(
    "gap (%d) is smaller than the %s ceiling (%d) -- this can change clustering results (silhouette/density) or over-constrain each chunk's search (dynamic). Set gap >= %d.",
    gap, if (clustering_method == "dynamic") "cut_height" else "h_max", .height_ceiling, .height_ceiling
  ))
}

# ---- renumber cluster ids 1..k by ascending minimum position --------------
# cutree()'s raw label order is not guaranteed to follow spatial order, and
# cutreeDynamic's is even less likely to (see header note). Shared by all
# three methods so every chunk's group ids increase monotonically with
# position, and downstream group_offset math in cluster_chromosome() stays
# exact (labels become contiguous 1..k as a side effect).
renumber_by_position <- function(positions, labels) {
  min_pos_by_cluster <- tapply(positions, labels, min)
  ordered_old_ids     <- names(sort(min_pos_by_cluster))
  new_id_of           <- setNames(seq_along(ordered_old_ids), ordered_old_ids)
  unname(new_id_of[as.character(labels)])
}

# ---- clustering for one gap-chunk of positions: silhouette method ---------
# Builds one dist()/hclust() over this chunk's unique positions and picks the
# cut height maximizing mean silhouette width. Returns chunk-LOCAL group ids.
cluster_one_chunk_silhouette <- function(maf, pos_col, dist_method, hclust_method,
                                         h_min, h_max) {
  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  n <- length(positions)

  if (n < 2) {
    maf$group <- if (n == 0) integer(0) else 1L
    return(list(maf = maf, best_h = NA_integer_, best_score = NA_real_, n_positions = n))
  }

  d  <- dist(positions, method = dist_method)
  hc <- hclust(d, method = hclust_method)

  h_values  <- seq.int(h_min, h_max)
  sil_means <- rep(NA_real_, length(h_values))

  for (i in seq_along(h_values)) {
    lab <- tryCatch(cutree(hc, h = h_values[i]), error = function(e) NULL)
    if (is.null(lab)) next
    if (length(unique(lab)) > 1 && length(unique(lab)) < length(lab)) {
      sil_means[i] <- tryCatch(
        mean(silhouette(lab, d)[, "sil_width"]),
        error = function(e) NA_real_
      )
    }
  }

  if (all(is.na(sil_means))) {
    best_h     <- h_min
    best_score <- NA_real_
  } else {
    best_idx   <- which.max(sil_means)
    best_h     <- h_values[best_idx]
    best_score <- sil_means[best_idx]
  }

  labels <- tryCatch(cutree(hc, h = best_h),
                     error = function(e) seq_along(positions))   # fallback: singletons
  labels <- renumber_by_position(positions, labels)

  lab_df <- setNames(data.frame(positions, labels), c(pos_col, "group"))
  maf <- left_join(maf, lab_df, by = pos_col)

  list(maf = maf, best_h = best_h, best_score = best_score, n_positions = n)
}

# ---- clustering for one gap-chunk of positions: density method ------------
# Picks the cut height maximizing a density/fragmentation-penalized score
# (see header). Math unchanged from Houman's original.
cluster_one_chunk_density <- function(maf, pos_col, dist_method, hclust_method,
                                      h_min, h_max) {
  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  n <- length(positions)

  if (n < 2) {
    maf$group <- if (n == 0) integer(0) else 1L
    return(list(maf = maf, best_h = NA_integer_, best_score = NA_real_, n_positions = n))
  }

  d           <- dist(positions, method = dist_method)
  hc          <- hclust(d, method = hclust_method)
  dist_matrix <- as.matrix(d)

  # Positions are fixed for this chunk, so the median inter-position gap
  # (part of the score's denominator) is constant across all candidate heights.
  med_gap <- median(diff(positions))

  h_values <- seq.int(h_min, h_max)
  scores   <- rep(NA_real_, length(h_values))

  for (i in seq_along(h_values)) {
    lab <- tryCatch(cutree(hc, h = h_values[i]), error = function(e) NULL)
    if (is.null(lab)) next

    clusters <- unique(lab)
    if (length(clusters) <= 1) next

    span <- tapply(positions, lab, function(p) max(p) - min(p) + 1L)

    cl_score <- numeric(length(clusters))
    names(cl_score) <- as.character(clusters)
    for (cl in clusters) {
      idx <- which(lab == cl)
      k   <- length(idx)
      if (k > 1) {
        sub <- dist_matrix[idx, idx]
        cl_score[as.character(cl)] <-
          (k / span[as.character(cl)]) *
          (1 / (1 + mean(sub[lower.tri(sub)])))
      } else {
        cl_score[as.character(cl)] <- 0  # singleton: no pairwise distance
      }
    }

    # Penalise fragmentation: mean per-cluster score divided by the cluster
    # count raised to log2(median gap). (lambda is fixed at 1.)
    scores[i] <- mean(cl_score) / length(clusters)^log2(med_gap)
  }

  if (all(is.na(scores))) {
    best_h     <- h_min
    best_score <- NA_real_
  } else {
    best_idx   <- which.max(scores)
    best_h     <- h_values[best_idx]
    best_score <- scores[best_idx]
  }

  labels <- tryCatch(cutree(hc, h = best_h),
                     error = function(e) seq_along(positions))   # fallback: singletons
  labels <- renumber_by_position(positions, labels)

  lab_df <- setNames(data.frame(positions, labels), c(pos_col, "group"))
  maf <- left_join(maf, lab_df, by = pos_col)

  list(maf = maf, best_h = best_h, best_score = best_score, n_positions = n)
}

# ---- clustering for one gap-chunk of positions: dynamic method ------------
# Lets cutreeDynamic() (method = "hybrid") resolve cluster boundaries
# locally, branch by branch, instead of picking one height for the chunk.
cluster_one_chunk_dynamic <- function(maf, pos_col, dist_method, hclust_method,
                                      cut_height, min_cluster_size, deep_split,
                                      max_abs_core_scatter, min_abs_gap) {
  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  n <- length(positions)

  if (n < 2) {
    maf$group <- if (n == 0) integer(0) else 1L
    return(list(maf = maf,
                height_df = data.frame(group = integer(0), eff_height = numeric(0)),
                n_positions = n))
  }

  d  <- dist(positions, method = dist_method)
  hc <- hclust(d, method = hclust_method)

  # maxAbsCoreScatter/minAbsGap are built up via a named list and passed
  # through do.call() rather than as direct arguments, so that leaving them
  # unset (NULL) truly OMITS the argument -- letting cutreeDynamic compute
  # its own internal default -- rather than passing an explicit NULL, which
  # a function keying off missing() rather than is.null() could handle
  # differently than true omission.
  cd_args <- list(
    dendro         = hc,
    distM          = as.matrix(d),
    method         = "hybrid",
    minClusterSize = min_cluster_size,
    deepSplit      = deep_split,
    cutHeight      = cut_height,   # decoupled from `gap` -- see header note
    verbose        = 0
  )
  if (!is.null(max_abs_core_scatter)) cd_args$maxAbsCoreScatter <- max_abs_core_scatter
  if (!is.null(min_abs_gap))          cd_args$minAbsGap         <- min_abs_gap

  labels <- tryCatch(
    do.call(cutreeDynamic, cd_args),
    error = function(e) seq_along(positions)   # fallback: singletons
  )

  # cutreeDynamic labels unassigned ("grey"/noise) points 0 rather than
  # leaving them in their own cluster. Give each of those its own singleton
  # focus instead of silently merging unrelated noise points under label 0.
  if (any(labels == 0)) {
    noise_idx <- which(labels == 0)
    labels[noise_idx] <- max(labels) + seq_along(noise_idx)
  }

  labels <- renumber_by_position(positions, labels)

  # Per-focus "effective height": since cutreeDynamic never picks one height
  # for the whole chunk, best_h has no equivalent here. Instead, for each
  # resulting cluster, the height of the smallest dendrogram subtree that
  # spans all of its members equals the maximum cophenetic distance among
  # those members -- a well-defined per-focus stand-in for "what height was
  # this focus effectively cut at".
  coph <- as.matrix(cophenetic(hc))
  clusters <- unique(labels)
  eff_height <- vapply(clusters, function(cl) {
    idx <- which(labels == cl)
    if (length(idx) > 1) max(coph[idx, idx]) else NA_real_
  }, numeric(1))
  height_df <- data.frame(group = clusters, eff_height = eff_height)

  lab_df <- setNames(data.frame(positions, labels), c(pos_col, "group"))
  maf <- left_join(maf, lab_df, by = pos_col)

  list(maf = maf, height_df = height_df, n_positions = n)
}

# ---- dispatch to the configured method -------------------------------------
cluster_one_chunk <- function(maf, pos_col, dist_method, hclust_method, clustering_method,
                              h_min, h_max, cut_height, min_cluster_size, deep_split,
                              max_abs_core_scatter, min_abs_gap) {
  switch(clustering_method,
    silhouette = cluster_one_chunk_silhouette(maf, pos_col, dist_method, hclust_method,
                                              h_min, h_max),
    density    = cluster_one_chunk_density(maf, pos_col, dist_method, hclust_method,
                                           h_min, h_max),
    dynamic    = cluster_one_chunk_dynamic(maf, pos_col, dist_method, hclust_method,
                                           cut_height, min_cluster_size, deep_split,
                                           max_abs_core_scatter, min_abs_gap)
  )
}

# ---- split unique positions into gap-delimited chunks ----------------------
# A chunk boundary is placed between two consecutive unique positions whenever
# they are more than `gap` apart. For silhouette/density, this is exact, not
# approximate: for any two point sets that lie entirely on either side of
# such a gap G, every pairwise distance between a point on the left and a
# point on the right is >= G (by construction of 1D sorted coordinates), so
# the merge height for single, complete, average or centroid linkage -- and
# the centroid-distance term ward.D2 is built from -- between any left-subset
# and any right-subset is also >= G. Since G > gap >= h_max, no cut at
# h <= h_max can ever join points from different chunks, so clustering each
# chunk independently is equivalent to clustering the whole chromosome at any
# h <= h_max. For the dynamic method this proof does NOT transfer -- see the
# header note. Returns a chunk id per element of `positions` (already
# assumed sorted, unique).
split_into_gap_chunks <- function(positions, gap) {
  n <- length(positions)
  if (n < 2) return(rep(1L, n))
  gaps <- diff(positions)
  cumsum(c(TRUE, gaps > gap))
}

# ---- cluster a whole chromosome by clustering each gap-chunk separately ----
cluster_chromosome <- function(maf, pos_col, dist_method, hclust_method, clustering_method,
                               gap, h_min, h_max, cut_height, min_cluster_size, deep_split,
                               max_abs_core_scatter, min_abs_gap) {
  empty_plot_df <- function() {
    if (clustering_method == "dynamic") {
      data.frame(group = integer(0), eff_height = numeric(0))
    } else {
      data.frame(chunk = integer(0), best_h = integer(0),
                best_score = numeric(0), n_positions = integer(0))
    }
  }

  if (nrow(maf) == 0) {
    maf$group <- integer(0)
    return(list(maf = maf, plot_df = empty_plot_df(), n_chunks = 0L))
  }

  maf <- maf %>% arrange(.data[[pos_col]])
  positions <- sort(unique(maf[[pos_col]]))
  chunk_of_position <- split_into_gap_chunks(positions, gap)
  # Numeric index match (not a character-keyed lookup): positions here are
  # exactly the values maf[[pos_col]] was derived from, so no risk of a
  # scientific-notation mismatch on large genomic coordinates.
  maf$.chunk <- chunk_of_position[match(maf[[pos_col]], positions)]

  n_chunks <- max(chunk_of_position)
  message(sprintf(
    "  %d unique position(s) split into %d gap-chunk(s) (gap = %d, clustering_method = %s); largest chunk = %d position(s)",
    length(positions), n_chunks, gap, clustering_method, max(table(chunk_of_position))
  ))

  maf_chunks   <- vector("list", n_chunks)
  diag_rows    <- vector("list", n_chunks)
  group_offset <- 0L

  for (cid in seq_len(n_chunks)) {
    sub <- maf[maf$.chunk == cid, setdiff(names(maf), ".chunk"), drop = FALSE]
    res <- cluster_one_chunk(sub, pos_col, dist_method, hclust_method, clustering_method,
                             h_min, h_max, cut_height, min_cluster_size, deep_split,
                             max_abs_core_scatter, min_abs_gap)

    local_max <- if (nrow(res$maf) > 0) max(res$maf$group) else 0L
    res$maf$group <- res$maf$group + group_offset
    maf_chunks[[cid]] <- res$maf

    if (clustering_method == "dynamic") {
      if (nrow(res$height_df) > 0) {
        diag_rows[[cid]] <- data.frame(group = res$height_df$group + group_offset,
                                       eff_height = res$height_df$eff_height)
      }
    } else {
      diag_rows[[cid]] <- data.frame(chunk = cid, best_h = res$best_h,
                                     best_score = res$best_score,
                                     n_positions = res$n_positions)
    }

    group_offset <- group_offset + local_max
  }

  plot_df <- bind_rows(diag_rows)
  if (nrow(plot_df) == 0) plot_df <- empty_plot_df()

  list(maf = bind_rows(maf_chunks), plot_df = plot_df, n_chunks = n_chunks)
}

# ---- run -------------------------------------------------------------------
maf_chr <- read_tsv(maf_path, show_col_types = FALSE, progress = FALSE)
message(sprintf("chr %s: %d rows, %d unique positions (clustering_method = %s)",
                chrom, nrow(maf_chr), length(unique(maf_chr[[pos_col]])), clustering_method))

res <- cluster_chromosome(maf_chr, pos_col, dist_method, hclust_method, clustering_method,
                          gap, h_min, h_max, cut_height, min_cluster_size, deep_split,
                          max_abs_core_scatter, min_abs_gap)

write_tsv(res$maf, out_tsv, na = "")

if (make_plots) {
  if (clustering_method == "dynamic") {
    # One point per FOCUS (not per chunk, not per height): cutreeDynamic has
    # no single best_h per chunk, so this shows each resulting focus's own
    # effective height instead.
    plot_data <- dplyr::filter(res$plot_df, !is.na(eff_height))
    p <- ggplot(plot_data, aes(x = eff_height)) +
      geom_histogram(binwidth = 1, boundary = -0.5, fill = "grey40", colour = NA) +
      geom_vline(xintercept = cut_height, linetype = "dashed", colour = "firebrick") +
      labs(x = "Effective cut height per focus (max cophenetic distance within it)",
           y = "Number of foci",
           title = sprintf("Distribution of effective cut height (dynamic) — chr %s (%d gap-chunk(s))",
                           chrom, res$n_chunks),
           subtitle = "Dashed line: cut_height ceiling. Singleton foci (no internal distance) are excluded.") +
      theme_minimal()
    ggsave(out_plot, p, width = 7, height = 4)

    # No native score for this method -- see header note.
    placeholder <- ggplot() + theme_void() +
      labs(title = sprintf("No native score for clustering_method = \"dynamic\" — chr %s", chrom),
           subtitle = "cutreeDynamic never computes a scalar goodness score; see cluster_foci.R header.")
    ggsave(out_score_plot, placeholder, width = 7, height = 4)

  } else {
    method_label <- if (clustering_method == "silhouette") {
      "Mean silhouette width"
    } else {
      "Density score (no fixed scale -- relative within this run only)"
    }

    # One point per gap-chunk: shows whether best_h is piling up against
    # h_max across the chromosome, the actionable signal for whether
    # h_min/h_max are sized well.
    plot_data <- dplyr::filter(res$plot_df, !is.na(best_h))
    p <- ggplot(plot_data, aes(x = best_h)) +
      geom_histogram(binwidth = 1, boundary = h_min - 0.5, fill = "grey40", colour = NA) +
      geom_vline(xintercept = h_min, linetype = "dashed", colour = "steelblue") +
      geom_vline(xintercept = h_max, linetype = "dashed", colour = "firebrick") +
      labs(x = "Chosen cut height (best_h) per gap-chunk",
           y = "Number of gap-chunks",
           title = sprintf("Distribution of chosen cut height (%s) — chr %s (%d gap-chunk(s))",
                           clustering_method, chrom, res$n_chunks),
           subtitle = "Dashed lines: h_min (blue) / h_max (red). Mass piling up at h_max suggests the search range is too narrow.") +
      theme_minimal()
    ggsave(out_plot, p, width = 7, height = 4)

    # One point per gap-chunk: its chosen height against the score it
    # achieved there.
    score_data <- dplyr::filter(res$plot_df, !is.na(best_h) & !is.na(best_score))
    p_score <- ggplot(score_data, aes(x = best_h, y = best_score)) +
      geom_point(alpha = 0.5, colour = "grey20") +
      geom_vline(xintercept = h_min, linetype = "dashed", colour = "steelblue") +
      geom_vline(xintercept = h_max, linetype = "dashed", colour = "firebrick") +
      labs(x = "Chosen cut height (best_h) per gap-chunk",
           y = method_label,
           title = sprintf("Cut-height confidence (%s) — chr %s (%d gap-chunk(s))",
                           clustering_method, chrom, res$n_chunks),
           subtitle = "Each point is one gap-chunk. Low scores mean the chosen height wasn't a strong/confident peak.") +
      theme_minimal()
    ggsave(out_score_plot, p_score, width = 7, height = 4)
  }
} else {
  placeholder <- ggplot() + theme_void() +
    labs(title = sprintf("Plotting disabled (make_plots: False) — chr %s", chrom))
  ggsave(out_plot, placeholder, width = 7, height = 4)
  ggsave(out_score_plot, placeholder, width = 7, height = 4)
}

message(sprintf("chr %s: %d gap-chunk(s), %d foci total -> %s",
                chrom, res$n_chunks, length(unique(res$maf$group)), out_tsv))
