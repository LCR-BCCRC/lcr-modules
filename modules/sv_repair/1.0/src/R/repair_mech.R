#!/usr/bin/env Rscript

## merge vcf dataframe with rmsk and non_b dna; identify repair mechanism

log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressWarnings(
  suppressPackageStartupMessages(
    library(tidyverse)
  )
)

# Read in the vcf_df
vcf <- suppressMessages(
  read_tsv(
    snakemake@input[[1]],
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/02-vcf_to_df/genome--hg38/11-19987T--11-19987N--matched/somaticSV.pass.precise.tsv",
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/02-vcf_to_df/genome--hg38/09-27079T--09-27079N--matched/somaticSV.pass.precise.tsv",
    comment="#",
    col_names=TRUE
  )
 %>%
  rename(end = pos)
) %>% suppressMessages() %>% suppressWarnings()

# Read in the rmsk
rmsk <- read_tsv(
    snakemake@input[[2]],
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/04-rmsk/genome--hg38/11-19987T--11-19987N--matched/intersect.rmsk.bed",
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/04-rmsk/genome--hg38/09-27079T--09-27079N--matched/intersect.rmsk.bed",
    comment="#",
    col_names = c("chr","start","end","chr.feature","start.feature","end.feature","name","length","strand","v1","v2","v3","v4","real.name", "v5","v6","v7","v8")
  ) %>%
  group_by(chr, start, end) %>%
  summarize(feature_rmsk = paste0(real.name, collapse=","), name_rmsk = paste0(name, collapse=",")) %>% 
  suppressMessages() %>% suppressWarnings()

# Read in the non_b
non_b <-  read_tsv(
    snakemake@input[[3]],
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/05-nonb_dna/genome--hg38/11-19987T--11-19987N--matched/intersect.nonb.bed",
#    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/05-nonb_dna/genome--hg38/09-27079T--09-27079N--matched/intersect.nonb.bed",
    comment="#",
    col_names = c("chr","start","end","chr.feature","start.feature","end.feature","name","length","strand","v1","feature.name","v3","deets")
) %>% 
  mutate(
    sequence = str_extract(deets, "sequence=[a-z]+")
  ) %>%
  group_by(chr, start, end) %>%
  summarize(feature_nonb = paste0(feature.name, collapse=","), nonb_seq = paste0(sequence, collapse=","))%>%
  separate_wider_delim(., cols=nonb_seq, delim=",", names_sep="", too_few = "align_start") %>% suppressMessages() %>% suppressWarnings()


#ifelse((nrow(non_b) == nrow(vcf) & nrow(rmsk) == nrow(vcf)), NA, "STOP")

merged.1 <- merge( vcf, rmsk, by.x=c("chrom","end"), by.y=c("chr","end"))
merged <- merge( merged.1, non_b, by.x=c("chrom","end"), by.y = c("chr","end")) %>% 
  filter(SVTYPE == "BND") %>%
  mutate_if(is.character, list(~na_if(., ".")))

##group by mates??
df_id_mate <- merged

if (nrow(merged)>0) {
  df_id_mate <- merged %>% 
    rowwise() %>% 
    mutate(ID_mate = str_c(min(ID, MATEID), max(ID, MATEID))) 
}

paired <- df_id_mate
if (nrow(merged) == 0) {
  summary.repair <- data.frame("x" = c(0,0,0,0,0,0))
  rownames(summary.repair)<- c("unknown","cNHEJ","possible MMEJ","probable MMEJ", "possible HR","probable HR")
}
  
  
if (nrow(merged)>0) {
paired <- df_id_mate %>% 
  left_join(
    df_id_mate, 
    by = c(
      "ID_mate", 
      "ID" = "MATEID", 
      "MATEID" = "ID"
    ), 
    suffix = c("_ID", "_mate")
  ) %>% 
  group_by(ID_mate) %>% 
  slice_min(ID, n=1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(-ID_mate) %>%
  mutate_if(is.character, list(~na_if(., "NA")))


x <- c( "nonb_seq1_ID","nonb_seq1_mate","nonb_seq2_mate","nonb_seq3_mate")

paired[x[!(x %in% colnames(paired))]] = 0


paired <- paired %>% mutate(
  repair_mech_hom = case_when(
    is.na(HOMLEN_ID) ~ "cNHEJ",
    HOMLEN_ID < 3 ~ "cNHEJ",
    HOMLEN_ID == 3 ~ "possible MMEJ",
    HOMLEN_ID >3 & HOMLEN_ID < 31 ~ "probable MMEJ",
    HOMLEN_ID > 30 ~ "probable HR",
    TRUE ~ "no repair"),
  repair_mech_feature = case_when(
    name_rmsk_ID == name_rmsk_mate ~ "possible HR",
    nonb_seq1_ID == nonb_seq1_mate~ "possible HR",
    nonb_seq1_ID == nonb_seq2_mate~ "possible HR",
    nonb_seq1_ID == nonb_seq3_mate~ "possible HR",
    TRUE ~ "no repair"),
  repair_summary = case_when(
    repair_mech_hom == "possible MMEJ" ~ "possible MMEJ",
    repair_mech_hom == "probable MMEJ" ~ "probable MMEJ",
    repair_mech_feature == "possible HR" ~ "possible HR",
    repair_mech_hom == "probable HR" ~ "probable HR",
    repair_mech_hom == "cNHEJ" ~ "cNHEJ",
    TRUE ~ "unknown"
  )
  
) 

summary.repair <- factor(paired$repair_summary, levels=c("unknown","cNHEJ","possible MMEJ","probable MMEJ","possible HR", "probable HR")) %>% summary() %>% as.data.frame()

}

colnames(summary.repair) <- snakemake@params[[1]]
write.table(summary.repair, snakemake@output[[1]],  sep="\t", quote = FALSE)
write.table(paired, snakemake@output[[2]], sep="\t", quote=FALSE, row.names = FALSE)

sink()
