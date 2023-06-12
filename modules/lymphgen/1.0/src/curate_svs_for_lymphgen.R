source("/projects/rmorin/projects/gambl-repos/gambl-lhilton/src/R/GAMBLR_libs.R")
library(GAMBLR)

# different cases for inputs
#   1) fish w/ sample_id + svar master -- has_fish & has_sv
#   2) fish w/ sample_id + manta -- has_fish & has_sv
#   3) fish w/o sample_id (like empty fish) + svar master -- ~has_fish & has_sv
#   4) fish w/o sample_id (like empty fish) + manta -- ~has_fish & has_sv
#   5) fish w/ sample_id + empty sv -- has_fish & ~has_sv

all_fish <- snakemake@input[["fish"]]
sv <- snakemake@input[["sv"]]

if (str_detect(sv, "manta")) { # manta
  sv <- fread(sv, skip = 130) %>% 
    as.data.frame() %>% 
    rename(CHROM_A = '#CHROM_A') %>% 
    select(CHROM_A, START_A, END_A, QUAL, STRAND_A, CHROM_B, START_B, END_B, STRAND_B) %>% 
    rename(SCORE = QUAL) %>% 
    mutate(tumour_sample_id = snakemake@wildcards[["tumour_id"]])
} else if (str_detect(sv, "svar_master")) { # svar_master
  sv <- fread(sv) %>% 
    as.data.frame() %>% 
    select(CHROM_A, START_A, END_A, SCORE, STRAND_A, CHROM_B, START_B, END_B, STRAND_B, tumour_sample_id)
} else { # empty
  sv <- fread(sv) %>% 
    as.data.frame()
}

all_fish <- fread(all_fish) %>% as.data.frame()
one_fish <- all_fish %>%
  filter(sample_id == snakemake@wildcards[["tumour_id"]])

config_bcl6 <- snakemake@config[["lcr-modules"]][["lymphgen"]][["options"]][["real_bcl6"]]
config_bcl2 <- snakemake@config[["lcr-modules"]][["lymphgen"]][["options"]][["real_bcl2"]]

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
  svar_master_annotated <- annotate_sv(sv) %>%
    filter(!str_detect(fusion, "^NA-")) %>%
    filter(str_detect(fusion, "BCL[26]"))
}

real_bcl6 <- svar_master_annotated %>%
  filter(str_detect(fusion, "BCL6")) %>%
  filter(str_detect(fusion, paste(config_bcl6, collapse = "|"))) %>%
  pull(fusion)

real_bcl2 <- svar_master_annotated %>%
  filter(str_detect(fusion, "BCL2")) %>%
  filter(str_detect(fusion, paste(config_bcl2, collapse = "|"))) %>%
  pull(fusion)

if (nrow(one_fish) == 0) {
  final_sv <- one_fish %>%
    add_row(sample_id = snakemake@wildcards[["tumour_id"]],
            BCL2_BA_consensus = case_when(!(is_empty(real_bcl2)) ~ "POS",
      TRUE ~ "NEG"), BCL6_BA_consensus = case_when(!(is_empty(real_bcl6)) ~ "POS",
        TRUE ~ "NEG"))
} else {
  final_sv <- one_fish %>%
    mutate(BCL2_BA_consensus = case_when(
      BCL2_BA_consensus == "POS" | !(is_empty(real_bcl2)) ~ "POS",
      TRUE ~ "NEG"
    )) %>%
    mutate(BCL6_BA_consensus = case_when(
      BCL6_BA_consensus == "POS" | !(is_empty(real_bcl6)) ~ "POS",
      TRUE ~ "NEG"
    ))
}

write_tsv(final_sv, snakemake@output[[1]])
