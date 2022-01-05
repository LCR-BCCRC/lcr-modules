suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
  })
  
  final_cols <- c(
    "CHROM_A", 
    "START_A", 
    "END_A",
    "CHROM_B", 
    "START_B", 
    "END_B",
    "manta_name", 
    "SCORE", 
    "STRAND_A",
    "STRAND_B",
    "tumour_sample_id",
    "normal_sample_id",
    "VAF",
    "DP",
    "ANNOTATION_A", 
    "DIST_TO_ANNOTATION_A", 
    "ANNOTATION_B", 
    "DIST_TO_ANNOTATION_B",
    "gridss_name",
    "delly_name"
  )
  
  # Load and tidy combined GRIDSS and Manta bedpe annotations
  combined_a <- read_tsv(snakemake@input[["a_combined"]], 
                         col_types = cols(
                           CHROM_A = col_character(),
                           CHROM_B = col_character(),
                           VAF = col_double()
                         ),
                         col_names <- c(
                           "CHROM_A", 
                           "START_A", 
                           "END_A", 
                           "manta_name",
                           "SCORE",
                           "STRAND_A",
                           "STRAND_B",
                           "tumour_sample_id",
                           "normal_sample_id",
                           "VAF",
                           "DP",
                           "gridss_name", 
                           "CHROM_ANN_A", 
                           "START_ANN_A", 
                           "END_ANN_A", 
                           "ANNOTATION_A", 
                           "DIST_TO_ANNOTATION_A"
                         ))
  
  combined_b <- read_tsv(snakemake@input[["b_combined"]], 
                         col_types = cols(
                           CHROM_A = col_character(),
                           CHROM_B = col_character(),
                           VAF = col_double()
                         ), 
                         col_names <- c(
                           "CHROM_B", 
                           "START_B", 
                           "END_B", 
                           "manta_name",
                           "SCORE",
                           "STRAND_A",
                           "STRAND_B",
                           "tumour_sample_id",
                           "normal_sample_id",
                           "VAF",
                           "DP",
                           'gridss_name', 
                           "CHROM_ANN_B", 
                           "START_ANN_B", 
                           "END_ANN_B", 
                           "ANNOTATION_B", 
                           "DIST_TO_ANNOTATION_B"
                         ))
  
  
  
  if(nrow(combined_a) > 0){
    combined_bedpe <- full_join(combined_a, combined_b)
    combined_bedpe <- combined_bedpe %>% 
      # Select useful columns
      select(c(
        "CHROM_A", 
        "START_A", 
        "END_A",
        "CHROM_B", 
        "START_B", 
        "END_B",
        "manta_name",
        "SCORE",
        "STRAND_A",
        "STRAND_B",
        "tumour_sample_id",
        "normal_sample_id",
        "VAF",
        "DP",
        "gridss_name", 
        "ANNOTATION_A", 
        "DIST_TO_ANNOTATION_A", 
        "ANNOTATION_B", 
        "DIST_TO_ANNOTATION_B"
      )) %>% 
      # Remove small intrachromosomsal events
      filter(!(CHROM_A == CHROM_B & abs(START_A - START_B) < 5000)) %>% 
      distinct(CHROM_A, CHROM_B, START_A, START_B, END_A, END_B, STRAND_A, STRAND_B, .keep_all = TRUE)
  } 
write_tsv(combined_bedpe, snakemake@output[["bedpe"]])

save.image(snakemake@log[[1]])

})