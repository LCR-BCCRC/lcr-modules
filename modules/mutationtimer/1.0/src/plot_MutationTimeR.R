#!/usr/bin/env Rscript

# Script to plot the timing info from MutationTimeR for SSMs and CNAs
# meant to be run in the mutationtimer lcr-module
# output file format should be pdf!!

# Usage:
# Rscript plot_MutationTimeR.R <path/to/mutationtimer/timed_ssm.tsv> <path/to/mutationtimer/timed_cna.tsv> <path/for/output/full_plot.png> <path/for/output/min_plot.png> <tumour_sample_id> <normal_sample_id> <projection> <path/for/log>


suppressWarnings(
suppressPackageStartupMessages({
  library(BSgenome)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GAMBLR.helpers)
  library(data.table)
  library(ggpubr)
  library(readr)
  library(stringr)
  library(dplyr)
  library(ggplot2)
})
)

select = dplyr::select
filter = dplyr::filter
count = dplyr::count

# Parse command-line arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("input_ssm", "input_cna", "output_full", "output_min", "tumour_sample_id", "normal_sample_id", "projection", "log_path")
args <- setNames(args, arg_names[1:length(args)])

# Log both the stdout and stderr
log <- file(args$log_path, open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Print args for de-bugging -----------------------------------------------------------
cat(paste("timed ssm file:", args$input_ssm, "\n"))
cat(paste("timed cna file:", args$input_cna, "\n"))
cat(paste("path for full plot:", args$output_full, "\n"))
cat(paste("path for min plot:", args$output_min, "\n"))
cat(paste("tumour_sample_id:", args$tumour_sample_id, "\n"))
cat(paste("normal_sample_id:", args$normal_sample_id, "\n"))
cat(paste("projection:", args$projection, "\n"))
cat(paste("path for log:", args$log_path, "\n"))

# Set command-line argument values -----------------------------------------------------
input_ssm <- args$input_ssm
input_cna <- args$input_cna
output_full <- args$output_full
output_min <- args$output_min
tumour_id <- args$tumour_sample_id
normal_id <- args$normal_sample_id
projection <- args$projection

# Read in timed data
# -----------------------------------------------------
cat("Reading in timed SSM data...\n")
timed_ssm <- read_tsv(input_ssm, show_col_types = FALSE, guess_max = 3000) %>%
  mutate(Chromosome = as.character(Chromosome)) # otherwise non-prefixed ones get assigned as numeric when X is not present
cat("Reading in timed CNA data...\n")
timed_cna <- read_tsv(input_cna, show_col_types = FALSE) %>%
  mutate(chr = as.character(chr))

# Check for CNAs with time info
# -----------------------------------------------------
cat("Checking for CNAs with time info...\n")
if(dim(timed_cna %>% filter(!is.na(time)))[1] == 0){
  cat("No timed CNAs found. Min plot will be empty.\n")
}else {
  cat("Time info for CNAs found.\n")
}

# Plotting Functions
# -----------------------------------------------------
plot_timed_SSM <- function(timed_ssm,
                    timed_cna,
                    all_ssm = FALSE,
                    genome_build = "hg38",
                    base_size = 12,
                    point_size = 0.5){
  if(genome_build %in% "grch37"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg19
    #need to remove chr prefix for grch37
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
    names(chrom_lengths) <- str_remove(names(chrom_lengths), "chr")
  }else if(genome_build == "hg38"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg38
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
  }else{
    stop("Unknown genome specified")
  }
  s <- names(chrom_lengths)
  timed_ssm <- timed_ssm %>%
    mutate(VAF = t_alt_count/t_depth,
            start = Start_Position,
            end = End_Position,
            chr = Chromosome)
  timed_ssm$Chromosome = factor(timed_ssm$Chromosome, levels=s)
  timed_ssm$cls <- factor(timed_ssm$CLS_time_label, levels = c("clonal [early]", "clonal [late]", "clonal [NA]", "subclonal"))

  if(all_ssm){
    p <- ggplot(data = timed_ssm, aes(x=Start_Position, y=VAF, color=cls)) +
      geom_point(alpha=0.7, size=point_size, show.legend=TRUE) +
      facet_wrap(~Chromosome, scales="free_x", nrow=1) +
      ylim(c(0,1)) +
      theme_Morons(base_size = base_size) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
      ylab("VAF") +
      scale_colour_manual(values = c("clonal [early]"="#C77CFF", "clonal [late]"="#7CAE00", "clonal [NA]"="#00BFC4",  "subclonal"="#F8766D"),
      drop = FALSE)
  }else{
    just_timed <- timed_cna %>%
      filter(!is.na(time)) %>%
      mutate(start=startpos, end=endpos)
    x <- data.table(timed_ssm)
    y <- data.table(just_timed)

    setkey(y, chr, start, end)

    kept_ssm <- foverlaps(x, y, type="within") %>%
      dplyr::filter(!is.na(startpos)) %>%
      as.data.frame()

    p <- ggplot(data = kept_ssm, aes(x=Start_Position, y=VAF, color=cls)) +
      geom_point(alpha=0.7, size=point_size, show.legend=TRUE) +
      facet_wrap(~Chromosome, scales="free_x", nrow=1) +
      ylim(c(0,1)) +
      theme_Morons(base_size = base_size) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
      ylab("VAF") +
      scale_colour_manual(values = c("clonal [early]"="#C77CFF", "clonal [late]"="#7CAE00", "clonal [NA]"="#00BFC4",  "subclonal"="#F8766D"),
      drop = FALSE)
  }

  return(p)
}


plot_timed_cna <- function(timed_cna,
                          genome_build = "hg38",
                          base_size = 12,
                          all_chrom = FALSE){
  if(genome_build %in% "grch37"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg19
    #need to remove chr prefix for grch37
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
    names(chrom_lengths) <- str_remove(names(chrom_lengths), "chr")
  }else if(genome_build == "hg38"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg38
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
  }else{
      stop("Unknown genome specified")
  }
  s <- names(chrom_lengths)
  l <- as.numeric(chrom_lengths)
  cn_bounds <- data.frame(chr = s, start=1, end=l)

  timed_cna$chr <- factor(timed_cna$chr, levels=s)
  cn_bounds$chr <- factor(cn_bounds$chr, levels=s)

  if(all_chrom){
    if( dim(timed_cna %>% filter(!is.na(time)))[1] == 0){
      p <- ggplot() +
        ggtitle("No CNAs were timed")
    }else{
      p <- ggplot(timed_cna) +
        geom_segment(aes(y=time, yend=time, x=startpos, xend=endpos), colour="red")   +
        geom_rect(aes(xmin=endpos, xmax=startpos, ymin=time.lo, ymax=time.up), alpha=0.2) +
        geom_point(data=cn_bounds, aes(x=start,y=1), colour="white") +
        geom_point(data=cn_bounds, aes(x=end,y=1), colour="white") +
        ylim(c(0,1)) +
        facet_wrap(~chr, scales="free_x" ,nrow=1) +
        theme_Morons(base_size = base_size) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
        ylab("Time")
    }
  }else{
    just_timed <- timed_cna %>% filter(!is.na(time))
    cn_bounds <- cn_bounds %>% filter(chr %in% just_timed$chr)

    p <- ggplot(just_timed) +
      geom_segment(aes(y=time, yend=time, x=startpos, xend=endpos), colour="red") +
      geom_rect(aes(xmin=endpos, xmax=startpos, ymin=time.lo, ymax=time.up), alpha=0.2) +
      geom_point(data=cn_bounds, aes(x=start, y=1), colour="white") +
      geom_point(data=cn_bounds, aes(x=end, y=1), colour="white") +
      ylim(c(0,1)) +
      facet_wrap(~chr, scales="free_x", nrow=1) +
      theme_Morons(base_size = base_size) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
      ylab("Time")
  }

  return(p)
}



plot_total_CN <- function(timed_cna,
                          genome_build = "hg38",
                          base_size = 12,
                          all_chrom = FALSE){
  if(genome_build %in% "grch37"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg19
    #need to remove chr prefix for grch37
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
    names(chrom_lengths) <- str_remove(names(chrom_lengths), "chr")
  }else if(genome_build == "hg38"){
    bs_genome <- BSgenome.Hsapiens.UCSC.hg38
    chrom_lengths <- seqlengths(bs_genome)[c(1:23)] # only canonical and X
  }else{
    stop("Unknown genome specified")
  }
  s <- names(chrom_lengths)
  l <- as.numeric(chrom_lengths)

  timed_cna$chr <- factor(timed_cna$chr, levels=s)

  timed_cna_total <- timed_cna %>%
      mutate(size = endpos-startpos+1) %>%
      mutate(centre = startpos+size/2) %>%
      mutate(CN = major_cn+0.1, CN2 = minor_cn-0.1)

  just_timed_chrom <- timed_cna_total %>%
      filter(!is.na(time)) %>%
      pull(chr)
  just_timed <- timed_cna_total %>%
      filter(chr %in% just_timed_chrom)

  if(all_chrom){
    p <- ggplot(timed_cna_total) +
      geom_tile(aes(x=centre, y=CN, width=size, height=0.2, alpha=clonal_frequency), fill="lightblue") +
      geom_tile(aes(x=centre, y=CN2, width=size, height=0.2, alpha=clonal_frequency), fill="orange") +
      facet_wrap(~chr, scales="free_x", nrow=1) +
      ylim(c(-1,4)) +
      theme_Morons(base_size = base_size) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), axis.title.y=element_text(margin = margin(l=16)))
  }else{
    p <- ggplot(just_timed) +
      geom_tile(aes(x=centre, y=CN, width=size, height=0.2, alpha=clonal_frequency), fill="lightblue") +
      geom_tile(aes(x=centre, y=CN2, width=size, height=0.2, alpha=clonal_frequency), fill="orange") +
      facet_wrap(~chr, scales="free_x", nrow=1) +
      ylim(c(-1,4)) +
      theme_Morons(base_size = base_size) +
      theme(axis.text.x=element_blank() ,axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), axis.title.y=element_text(margin = margin(l=16)))
  }

    return(p)
}


plot_all <- function(timed_ssm,
                    timed_cna,
                    genome_build = "hg38",
                    sample_ids = "",
                    title = "",
                    base_size = 12,
                    point_size = 0.8,
                    all_ssm = FALSE,
                    all_chrom = FALSE){

  if(all_ssm & all_chrom){
    CN_total <- plot_total_CN(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom=TRUE)
    CN_timing <- plot_timed_cna(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom=TRUE)
    SSM_timing <- plot_timed_SSM(timed_ssm, timed_cna, all_ssm = TRUE, point_size = point_size, genome_build = genome_build, base_size = base_size)
    p <- ggarrange(SSM_timing, CN_timing, CN_total, ncol=1) %>%
      annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
  }else if(all_ssm != all_chrom){
    stop("all_ssm and all_chrom are not the same. Set both to either FALSE or TRUE.")
  }else if(dim(timed_cna %>% filter(!is.na(time)))[1] == 0){
    title <- "No CNAs were timed"
    p <- ggarrange(ggplot()) %>%
      annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
  }else{
    CN_total <- plot_total_CN(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom = FALSE)
    CN_timing <- plot_timed_cna(timed_cna, genome_build = genome_build, base_size = base_size, all_chrom = FALSE)
    SSM_timing <- plot_timed_SSM(timed_ssm, timed_cna, all_ssm = FALSE, point_size = point_size, genome_build = genome_build, base_size = base_size)
    p <- ggarrange(SSM_timing, CN_timing, CN_total, ncol=1) %>%
        annotate_figure(., top = text_grob(title, face = "bold", size = 12), fig.lab = sample_ids, fig.lab.pos = "top.left")
  }

  return(p)
}


# Creating and saving the plotting functions
# -----------------------------------------------------
cat("Making full plot...\n")
plot_full <- plot_all(timed_ssm, timed_cna, all_ssm = TRUE, all_chrom = TRUE, genome_build = projection,
  title = "All SSM and Copy Number States", sample_ids = paste0(tumour_id, "--", normal_id))

ggsave(filename = output_full, plot = plot_full, device = cairo_pdf, width = 12, height = 12, create.dir = TRUE)

cat("Making minimum plot...\n")
plot_min <- plot_all(timed_ssm, timed_cna, all_ssm = FALSE, all_chrom = FALSE, genome_build = projection,
  title ="SSM and Copy Number in Timed CNA Regions", sample_ids = paste0(tumour_id, "--", normal_id))

ggsave(filename = output_min, plot = plot_min, device = cairo_pdf, width = 12, height = 12, create.dir = TRUE)

cat("DONE!")
sink()