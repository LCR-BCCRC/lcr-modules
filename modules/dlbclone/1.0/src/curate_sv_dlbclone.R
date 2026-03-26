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
  make_option(c("-c","--real_bcl2"), type="character", help="Vector of real BCL2 fusion partners")
)

opt <- parse_args(OptionParser(option_list=option_list))

sv <- opt$sv_input
metadata <- read_tsv(opt$metadata)
output_metadata <- opt$output_metadata
config_bcl6 <- eval(parse(text = opt$real_bcl6))
config_bcl2 <- eval(parse(text = opt$real_bcl2))
config_myc <- c(".+")

# Read in data ------------------------------------------------------------

if (str_detect(Sys.readlink(sv), "empty")) { # empty
  print("empty")
  sv <- fread(sv) %>%
    as.data.frame()
} else if (str_detect(Sys.readlink(sv), "svar_master")) { # svar_master
  print("svar_master")
  sv <- fread(sv) %>%
    as.data.frame()# %>%
    # select(CHROM_A, START_A, END_A, SCORE, STRAND_A, CHROM_B, START_B, END_B, STRAND_B, tumour_sample_id)
} else { # manta
  print("not svar_master")
  sv <- fread(sv, skip = "CHROM") %>%
    #rename(CHROM_A = "#CHROM_A") %>%
    select(CHROM_A, START_A, END_A, SCORE, STRAND_A, CHROM_B, START_B, END_B, STRAND_B) %>%
    #rename(SCORE = QUAL) %>%
    mutate(tumour_sample_id = metadata$sample_id) 
}

sv$CHROM_A <- as.character(sv$CHROM_A)
sv$CHROM_B <- as.character(sv$CHROM_B)

if (nrow(sv) == 0) {
  svar_master_annotated <- data.frame(matrix(ncol = 14, nrow = 0))
  x <- c(
    "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score",
    "strand1", "strand2", "tumour_sample_id", "gene", "partner", "fusion"
  )
  colnames(svar_master_annotated) <- x
} else {
  svar_master_annotated <- annotate_sv(sv)
}

real_bcl6 <- svar_master_annotated %>%
  filter(str_detect(fusion, "BCL6")) %>%
  filter(str_detect(fusion, paste(config_bcl6, collapse = "|"))) %>%
  pull(fusion)

# print(real_bcl6)

real_bcl2 <- svar_master_annotated %>%
  filter(str_detect(fusion, "BCL2")) %>%
  filter(str_detect(fusion, paste(config_bcl2, collapse = "|"))) %>%
  pull(fusion)

# print(real_bcl2)

real_myc <- svar_master_annotated %>%
  filter(str_detect(fusion, "MYC")) %>%
  filter(!is.na(fusion)) %>%  # prevents NA from passing through
  filter(str_detect(fusion, paste(config_myc, collapse = "|"))) %>%
  pull(fusion)

# print(real_myc)

sample_id <- unique(sv$tumour_sample_id)[1]

BCL2_pos <- any(svar_master_annotated$partner %in% real_bcl2 & svar_master_annotated$gene == "BCL2")
BCL6_pos <- any(svar_master_annotated$partner %in% real_bcl6 & svar_master_annotated$gene == "BCL6")
MYC_pos <- any(svar_master_annotated$partner %in% real_myc & svar_master_annotated$gene == "MYC")

final_sv <- data.frame(
  sample_id = sample_id,
  BCL2_BA_consensus = NA_character_,
  BCL6_BA_consensus = NA_character_,
  MYC_BA_consensus = NA_character_
  ) %>%
    mutate(
      BCL2_BA_consensus = if_else(
        sample_id == !!sample_id,
        if_else(BCL2_pos, "POS", "NEG"),
        BCL2_BA_consensus
      ),
      BCL6_BA_consensus = if_else(
        sample_id == !!sample_id,
        if_else(BCL6_pos, "POS", "NEG"),
        BCL6_BA_consensus
      ),
      MYC_BA_consensus = if_else(
        sample_id == !!sample_id,
        if_else(MYC_pos, "POS", "NEG"),
        MYC_BA_consensus
      )
    )

# metadata combine SV
metadata <- metadata %>%
  inner_join(final_sv, by = "sample_id") 

write_tsv(metadata, output_metadata)
