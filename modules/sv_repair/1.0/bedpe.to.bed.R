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
    #snakemake@input[[1]],
    "/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair/02-window_bedpe/genome--hg38/01-20985T--01-20985N--matched/somaticSV.bedpe",
    comment="#",
    col_names=TRUE
  )
)

bedpe <- as.data.frame(
  bedpe
)
bedpe.A <- bedpe %>% select(X1, start1,end1,X7,X9)
bedpe.B <- bedpe %>% select(X4, start2,end2,X7,X10)
colnames(bedpe.A) <- c("chr","start","end","name","strand")
colnames(bedpe.B) <- c("chr","start","end","name","strand")

bed.file<-rbind(bedpe.A,bedpe.B)

# Saving the pre-processed data
write_tsv(bed.file, snakemake@output[[1]])

sink()
