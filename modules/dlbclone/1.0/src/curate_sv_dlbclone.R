#!/usr/bin/env Rscript
#
# Notes: 
#   This script is intended for use with the dlbclone-1.0 module in LCR-modules. 
# 
#   different cases for inputs
#     1) fish w/ sample_id + svar master -- has_fish & has_sv
#     2) fish w/ sample_id + manta -- has_fish & has_sv
#     3) fish w/o sample_id (like empty fish) + svar master -- ~has_fish & has_sv
#     4) fish w/o sample_id (like empty fish) + manta -- ~has_fish & has_sv
#     5) fish w/ sample_id + empty sv -- has_fish & ~has_sv

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
  library(GAMBLR.utils)
  library(optparse)
  library(tibble)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(readr)
})

# Determine arguments -----------------------------------------------------

# Define command line options
option_list <- list(
  make_option(c("-s", "--sv_input"), type = "character", help = "path to test sv file"),
  make_option(c("-m", "--metadata"), type = "character", help = "file path to metadata file"),
  make_option(c("-o", "--output_metadata"), type = "character", help = "output fusions file path"),
  make_option(c("-b","--real_bcl6"), type="character", help="Vector of real BCL6 fusion partners"),
  make_option(c("-c","--real_bcl2"), type="character", help="Vector of real BCL2 fusion partners"),
  make_option(c("-f","--sv_fish_colname"), type = "character", help = "A named vector that specifies the columns containing the FISH translocation status for any SV that is annotated in the metadata")
)

opt <- parse_args(OptionParser(option_list=option_list))

sv <- opt$sv_input
metadata <- read_tsv(opt$metadata)
output_metadata <- opt$output_metadata
config_bcl6 <- eval(parse(text = opt$real_bcl6))
config_bcl2 <- eval(parse(text = opt$real_bcl2))
config_myc <- c(".+")

parse_nullable <- function(x) {
    if (x == "None") {
        return(NULL)
    }
    eval(parse(text = x))
}

sv_fish_colname <- parse_nullable(opt$sv_fish_colname)

# Read in data ------------------------------------------------------------

if (str_detect(Sys.readlink(sv), "empty")) { # empty
  print("empty")
  sv <- fread(sv) %>%
    as.data.frame()
} else if (str_detect(Sys.readlink(sv), "svar_master")) { # svar_master
  print("svar_master")
  sv <- fread(sv) %>%
    as.data.frame()
} else { # manta
  print("not svar_master")
  sv <- fread(sv, skip = "CHROM")
  if (!"tumour_sample_id" %in% colnames(sv)) {
    stop("Manta SV missing tumour_sample_id column.")
  }
  sv <- sv %>%
    select(CHROM_A, START_A, END_A, SCORE, STRAND_A, CHROM_B, START_B, END_B, STRAND_B, tumour_sample_id)
}

sv$CHROM_A <- as.character(sv$CHROM_A)
sv$CHROM_B <- as.character(sv$CHROM_B)

# split per sample
sv_split <- split(sv, sv$tumour_sample_id)

results_list <- lapply(names(sv_split), function(sample_id) {

  sv_sub <- sv_split[[sample_id]]

  if (nrow(sv_sub) == 0) {
    svar_master_annotated <- data.frame(matrix(ncol = 14, nrow = 0))
    x <- c(
      "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score",
      "strand1", "strand2", "tumour_sample_id", "gene", "partner", "fusion"
    )
    colnames(svar_master_annotated) <- x
  } else {
    svar_master_annotated <- annotate_sv(sv_sub)
  }

  real_bcl6 <- svar_master_annotated %>%
    filter(str_detect(fusion, "BCL6")) %>%
    filter(str_detect(fusion, paste(config_bcl6, collapse = "|"))) %>%
    pull(fusion)

  real_bcl2 <- svar_master_annotated %>%
    filter(str_detect(fusion, "BCL2")) %>%
    filter(str_detect(fusion, paste(config_bcl2, collapse = "|"))) %>%
    pull(fusion)

  real_myc <- svar_master_annotated %>%
    filter(str_detect(fusion, "MYC")) %>%
    filter(!is.na(fusion)) %>%
    pull(fusion)

  BCL2_pos <- any(svar_master_annotated$partner %in% real_bcl2 & svar_master_annotated$gene == "BCL2")
  BCL6_pos <- any(svar_master_annotated$partner %in% real_bcl6 & svar_master_annotated$gene == "BCL6")
  MYC_pos  <- any(svar_master_annotated$partner %in% real_myc  & svar_master_annotated$gene == "MYC")

  data.frame(
    sample_id = sample_id,
    BCL2_BA_consensus = ifelse(BCL2_pos, "POS", "NEG"),
    BCL6_BA_consensus = ifelse(BCL6_pos, "POS", "NEG"),
    MYC_BA_consensus  = ifelse(MYC_pos,  "POS", "NEG")
  )
})

final_sv <- bind_rows(results_list)

# metadata combine SV
metadata <- metadata %>%
  inner_join(final_sv, by = "sample_id") 

# If FISH data is available for any SV, override the SV consensus call with the FISH call for any sample that is FISH POS. 
if (!is.null(sv_fish_colname)) {
  for (fish_col in sv_fish_colname) {

    gene <- sub("_F$", "", fish_col)
    ba_col <- paste0(gene, "_BA_consensus")

    if (ba_col %in% colnames(metadata) && fish_col %in% colnames(metadata)) {

      metadata[[ba_col]] <- ifelse(
        metadata[[fish_col]] == "POS",
        "POS",
        metadata[[ba_col]]
      )

    }
  }
}

write_tsv(metadata, output_metadata)
