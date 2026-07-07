#!/usr/bin/env Rscript

# Plot CCF output from _cnaqc_ccf rule
# Called via Snakemake script: directive; parameters accessed via snakemake@
# Inputs:  CCF TSV produced by run_ccf.R
# Outputs: PDF with per-chromosome VAF overview coloured by total CN

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(data.table)
  library(GAMBLR.helpers)
})

ccf_path  <- snakemake@input[["ccf"]]
out_plot  <- snakemake@output[["plot"]]
tumour_id <- snakemake@wildcards[["tumour_id"]]
width     <- snakemake@params[["width"]]
height    <- snakemake@params[["height"]]
base_size <- snakemake@params[["base_size"]]

cnaqc_overview_plot <- function(ccf_df, plot_type = "VAF", base_size = 6) {
  if (any(grepl("chr", ccf_df$chr))) {
    ccf_df <- mutate(ccf_df, chr = gsub("chr", "", chr))
  }
  ccf_df <- mutate(ccf_df,
    chr = ifelse(chr == "X", "23", chr),
    chr = as.numeric(chr)
  )

  cn_palette <- c(
    "2"  = "#2166AC",
    "1"  = "#74ADD1",
    "0"  = "#D1E5F0",
    "3"  = "#F4A582",
    "4"  = "#D6604D",
    "5"  = "#B2182B",
    "6+" = "#67001F"
  )

  cn_to_color <- function(cn) {
    key <- dplyr::case_when(
      is.na(cn) ~ NA_character_,
      cn >= 6   ~ "6+",
      TRUE      ~ as.character(cn)
    )
    unname(cn_palette[key])
  }

  seg_ranges <- ccf_df %>%
    filter(!is.na(bt_major_1)) %>%
    mutate(
      .major = bt_major_1,
      .minor = bt_minor_1,
      x_min  = as.integer(gsub("^[^:]+:([^:]+):.*$",       "\\1", segment_id)),
      x_max  = as.integer(gsub("^[^:]+:[^:]+:([^:]+):.*$", "\\1", segment_id))
    ) %>%
    distinct(chr, segment_id, .major, .minor, x_min, x_max) %>%
    tidyr::pivot_longer(c(.major, .minor), names_to = "allele", values_to = "cn") %>%
    mutate(
      allele    = if_else(allele == ".major", "Major", "minor"),
      seg_color = cn_to_color(cn),
      seg_lw    = if_else(allele == "Major", 4, 2)
    )

  y_col   <- if (plot_type == "depth") "DP" else "VAF"
  y_range <- range(ccf_df[[y_col]], na.rm = TRUE)
  y_span  <- diff(y_range)
  y_major <- y_range[2] + 0.08 * y_span
  y_minor <- y_range[2] + 0.03 * y_span

  seg_ranges <- mutate(seg_ranges,
    y_pos = if_else(allele == "Major", y_major, y_minor)
  )

  ccf_df <- ccf_df %>%
    mutate(
      total    = bt_major_1 + bt_minor_1,
      pt_color = cn_to_color(total),
      Position = from
    )

  all_cn      <- sort(unique(c(seg_ranges$cn[!is.na(seg_ranges$cn)],
                               ccf_df$total[!is.na(ccf_df$total)])))
  all_keys    <- ifelse(all_cn >= 6, "6+", as.character(all_cn))
  dedup       <- !duplicated(all_keys)
  legend_cols <- cn_to_color(all_cn[dedup])
  legend_labs <- all_keys[dedup]

  ggplot(ccf_df, aes(x = Position, y = .data[[y_col]],
                     shape = factor(mutation_multiplicity))) +
    geom_point(alpha = 0.7, aes(colour = pt_color)) +
    geom_segment(
      data = seg_ranges,
      aes(x = x_min, xend = x_max, y = y_pos, yend = y_pos,
          colour = seg_color, linewidth = seg_lw),
      inherit.aes = FALSE
    ) +
    scale_color_identity(
      name   = "Copy number",
      guide  = guide_legend(override.aes = list(linewidth = 6)),
      breaks = legend_cols,
      labels = legend_labs
    ) +
    scale_linewidth_identity(guide = "none") +
    scale_shape_discrete(name = "Multiplicity", na.translate = FALSE) +
    coord_cartesian(
      ylim = c(y_range[1], y_major + 0.02 * y_span),
      clip = "off"
    ) +
    facet_wrap(~chr, scales = "free_x") +
    theme_Morons(base_size = base_size)
}

# ---- Load and plot ----

ccf_df <- fread(ccf_path, sep = "\t", data.table = FALSE) %>% 
  mutate(mutation_multiplicity = ifelse(is.na(mutation_multiplicity), 0, mutation_multiplicity))
message(sprintf("Loaded %d mutations from %s", nrow(ccf_df), ccf_path))

pdf(out_plot, width = width, height = height)

if (all(is.na(ccf_df$bt_major_1))) {
  message("No CN-mapped mutations — writing blank plot page.")
  grid::grid.newpage()
  grid::grid.text(
    paste0(tumour_id, "\nNo CN-mapped mutations available for plotting"),
    gp = grid::gpar(fontsize = 14)
  )
} else {
  p <- cnaqc_overview_plot(ccf_df, base_size = base_size)
  p <- p + ggtitle(tumour_id)
  print(p)
}

dev.off()
message("Done.")
