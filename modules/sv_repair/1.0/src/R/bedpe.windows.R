#!/usr/bin/env Rscript

## manipulate bedpe file based on structural variant orientation to window of specified size

log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressWarnings(
  suppressPackageStartupMessages(
    library(tidyverse)
  )
)

# Read in the bedpe
bedpe <- suppressMessages(
  read_tsv(
    snakemake@input[[1]],
    #"/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair/01-input_bedpe/genome--hg38/01-20985T--01-20985N--matched/somaticSV.bedpe",
    comment="#",
    col_names=FALSE
   )
)
bedpe <- as.data.frame(
  bedpe
)

bedpe$strands <- paste0(bedpe$X9,bedpe$X10)

window.size<-snakemake@params[[1]]
#window.size<-1000
bedpe <- bedpe %>% 
  mutate(
    start1 = case_when(
      strands == "++" ~ X3-window.size,
      strands == "+-" ~ X3-window.size,
      strands == "-+" ~ X2-window.size,
      strands == "--" ~ X2-window.size
    ),
    end1 = case_when(
      strands == "++" ~ X3+window.size,
      strands == "+-" ~ X3+window.size,
      strands == "-+" ~ X2+window.size,
      strands == "--" ~ X2+window.size
    ),
    start2 = case_when(
      strands == "++" ~ X5-window.size,
      strands == "+-" ~ X6-window.size,
      strands == "-+" ~ X5-window.size,
      strands == "--" ~ X6-window.size
    ),
    end2 = case_when(
      strands == "++" ~ X5+window.size,
      strands == "+-" ~ X6+window.size,
      strands == "-+" ~ X5+window.size,
      strands == "--" ~ X6+window.size
    ),
  )

# Saving the pre-processed data
write_tsv(bedpe, snakemake@output[[1]])

sink()
