#!/usr/bin/env Rscript

# Genome-wide plot of one tumour's ModelSegments result: denoised copy ratio
# per target and alternate-allele fraction at heterozygous sites, each with
# the final segments. Targets removed before segmentation (outside the
# --min_log2/--max_log2 limits) are drawn in grey at the plot edges.
# Inputs:
#   - denoisedCR.tsv (unfiltered), hets.tsv, modelFinal.seg
# Outputs:
#   - PNG (or PDF, from the --out extension)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--denoised", type = "character", help = "denoisedCR.tsv from DenoiseReadCounts"),
  make_option("--hets", type = "character", help = "hets.tsv from ModelSegments"),
  make_option("--segments", type = "character", help = "modelFinal.seg from ModelSegments"),
  make_option("--sample", type = "character", help = "Title for the plot"),
  make_option("--out", type = "character", help = "Output file (.png or .pdf)"),
  make_option("--min_log2", type = "double", default = -5,
              help = "Lower limit used before segmentation [default: %default]"),
  make_option("--max_log2", type = "double", default = 5,
              help = "Upper limit used before segmentation [default: %default]"),
  make_option("--chromosome_regex", type = "character", default = "^(chr)?([1-9]|1[0-9]|2[0-2])$",
              help = "Contigs to plot [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
for (arg in c("denoised", "hets", "segments", "out")) {
  if (is.null(opt[[arg]])) stop(paste("Required argument missing:", arg))
}
if (is.null(opt$sample)) opt$sample <- sub("\\.modelFinal\\.seg$", "", basename(opt$segments))

read_gatk <- function(f) fread(cmd = paste("grep -v '^@'", shQuote(f)))

cr <- read_gatk(opt$denoised)
hets <- read_gatk(opt$hets)
seg <- read_gatk(opt$segments)

keep <- function(d) d[grepl(opt$chromosome_regex, as.character(CONTIG))]
cr <- keep(cr); hets <- keep(hets); seg <- keep(seg)

# cumulative genome coordinates, in natural chromosome order
chroms <- unique(as.character(cr$CONTIG))
chroms <- chroms[order(suppressWarnings(as.numeric(sub("^chr", "", chroms))), chroms)]
len <- cr[, .(len = max(END)), by = .(CONTIG = as.character(CONTIG))][match(chroms, CONTIG)]
len[, offset := c(0, cumsum(as.numeric(len))[-.N])]
off <- setNames(len$offset, len$CONTIG)
gpos <- function(chr, pos) off[as.character(chr)] + pos

# copy ratio: clip to the plotted range; removed targets in grey
ylim_cr <- c(-3, 3)
cr[, x := gpos(CONTIG, (START + END) / 2)]
cr[, removed := LOG2_COPY_RATIO < opt$min_log2 | LOG2_COPY_RATIO > opt$max_log2]
cr[, y := pmin(pmax(LOG2_COPY_RATIO, ylim_cr[1]), ylim_cr[2])]

# allele fraction of the tumour at the heterozygous sites
hets[, x := gpos(CONTIG, POSITION)]
hets[, y := ALT_COUNT / pmax(REF_COUNT + ALT_COUNT, 1)]

seg[, `:=`(x0 = gpos(CONTIG, START), x1 = gpos(CONTIG, END))]
seg_cr <- seg[, .(x0, x1, y = pmin(pmax(LOG2_COPY_RATIO_POSTERIOR_50, ylim_cr[1]), ylim_cr[2]),
                  track = "log2 copy ratio")]
# minor-allele fraction segments, mirrored at 1 - f to sit on both allele bands
maf <- seg[is.finite(MINOR_ALLELE_FRACTION_POSTERIOR_50)]
seg_af <- rbind(maf[, .(x0, x1, y = MINOR_ALLELE_FRACTION_POSTERIOR_50)],
                maf[, .(x0, x1, y = 1 - MINOR_ALLELE_FRACTION_POSTERIOR_50)])[, track := "allele fraction"]

points <- rbind(
  cr[, .(x, y, removed, track = "log2 copy ratio")],
  hets[, .(x, y, removed = FALSE, track = "allele fraction")]
)
tracks <- c("log2 copy ratio", "allele fraction")
points[, track := factor(track, tracks)]
seg_cr[, track := factor(track, tracks)]
seg_af[, track := factor(track, tracks)]

bounds <- data.table(x = c(len$offset[-1]))
mids <- len[, .(x = offset + len / 2, label = sub("^chr", "", CONTIG))]

p <- ggplot() +
  geom_vline(data = bounds, aes(xintercept = x), colour = "grey85", linewidth = 0.3) +
  geom_point(data = points[removed == FALSE], aes(x, y), size = 0.15, alpha = 0.25, colour = "grey30") +
  geom_point(data = points[removed == TRUE], aes(x, y), size = 0.3, alpha = 0.6, colour = "grey70") +
  geom_segment(data = seg_cr, aes(x = x0, xend = x1, y = y, yend = y), colour = "#C0392B", linewidth = 0.9) +
  geom_segment(data = seg_af[, .(x0, x1, y, track)], aes(x = x0, xend = x1, y = y, yend = y),
               colour = "#2E6FB7", linewidth = 0.9) +
  facet_grid(track ~ ., scales = "free_y") +
  scale_x_continuous(breaks = mids$x, labels = mids$label, expand = c(0.005, 0)) +
  labs(x = NULL, y = NULL,
       title = paste0(opt$sample, ": ", nrow(seg), " segments"),
       subtitle = sprintf("%d of %d targets removed before segmentation (grey, log2 outside %g to %g)",
                          sum(cr$removed), nrow(cr), opt$min_log2, opt$max_log2)) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95"))

ggsave(opt$out, p, width = 14, height = 6, dpi = 150)
