source("/projects/rmorin/projects/gambl-repos/gambl-lhilton/src/R/GAMBLR_libs.R")
library(GAMBLR)

all_fish <- snakemake@input[["fish"]]
one_svar_master <- snakemake@input[["svar_master"]]

one_svar_master <- fread(one_svar_master) %>% as.data.frame()
all_fish <- fread(all_fish) %>% as.data.frame()
one_fish <- all_fish %>%
  filter(sample_id == snakemake@wildcards[[3]])

config_bcl6 <- snakemake@config[["lcr-modules"]][["lymphgen"]][["options"]][["real_bcl6"]]
config_bcl2 <- snakemake@config[["lcr-modules"]][["lymphgen"]][["options"]][["real_bcl2"]]

one_svar_master$CHROM_A <- as.character(one_svar_master$CHROM_A)
one_svar_master$CHROM_B <- as.character(one_svar_master$CHROM_B)

if (nrow(one_svar_master) == 0) {
  svar_master_annotated <- data.frame(matrix(ncol = 14, nrow = 0))
  x <- c(
    "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score",
    "strand1", "strand2", "tumour_sample_id", "gene", "partner", "fusion"
  )
  colnames(svar_master_annotated) <- x
} else {
  one_svar_master_bcl <- one_svar_master %>%
    filter(str_detect(ANNOTATION_A, "BCL[26]") | str_detect(ANNOTATION_B, "BCL[26]"))

  svar_master_annotated <- annotate_sv(one_svar_master_bcl) %>%
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
  one_fish <- one_fish %>%
    mutate(across(everything(), as.character)) %>%
    add_row(sample_id = snakemake@wildcards[[3]], BCL2_BA_consensus = "NEG", BCL6_BA_consensus = "NEG")
}

final_sv <- one_fish %>%
  mutate(BCL2_BA_consensus = case_when(
    BCL2_BA_consensus == "POS" | !(is_empty(real_bcl2)) ~ "POS",
    TRUE ~ "NEG"
  )) %>%
  mutate(BCL6_BA_consensus = case_when(
    BCL6_BA_consensus == "POS" | !(is_empty(real_bcl6)) ~ "POS",
    TRUE ~ "NEG"
  ))

write_tsv(final_sv, snakemake@output[[1]])
