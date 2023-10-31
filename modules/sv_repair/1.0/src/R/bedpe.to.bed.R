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
    #"/projects/rmorin/projects/gambl-repos/gambl-kcoyle/results/gambl/sv_repair-1.0/01-window_bedpe/genome--grch37/14-11247T--14-11247N--matched/somaticSV.bedpe",
    comment="#",
    col_names=TRUE
  )
)

bedpe <- as.data.frame(
  bedpe
) %>% filter(X12 == "PASS")



bedpe.A <- bedpe %>% select(X1, start1,end1,X7,X9)
bedpe.B <- bedpe %>% select(X4, start2,end2,X7,X10)
colnames(bedpe.A) <- c("chr","start","end","name","strand")
colnames(bedpe.B) <- c("chr","start","end","name","strand")

bedpe.A<-bedpe.A %>% dplyr::group_by(name) %>% group_split()
bedpe.B<-bedpe.B %>% dplyr::group_by(name) %>% group_split()

write.tsv.funct.A <- function(x){
  obj.name <- x$name %>% unique()
  dir.name <- snakemake@output[[1]]
  dir.create(dir.name)
  file.name <- paste0(dir.name,"/",obj.name,".A.tsv")
  write_tsv(x, file = file.name, col_names=FALSE)
}

write.tsv.funct.B <- function(x){
  obj.name <- x$name %>% unique()
  dir.name <- snakemake@output[[1]]
  dir.create(dir.name)
  file.name <- paste0(dir.name,"/",obj.name,".B.tsv")
  write_tsv(x, file = file.name, col_names=FALSE)
}

bedpe.C <- bedpe %>% select(X1, X2,X3,X7,X9)
bedpe.D <- bedpe %>% select(X4, X5,X6,X7,X10)
colnames(bedpe.C) <- c("chr","start","end","name","strand")
colnames(bedpe.D) <- c("chr","start","end","name","strand")

bedpe.C<-bedpe.C %>% dplyr::group_by(name) %>% group_split()
bedpe.D<-bedpe.D %>% dplyr::group_by(name) %>% group_split()


write.tsv.funct.C <- function(x){
  obj.name <- x$name %>% unique()
  dir.name <- snakemake@output[[1]]
  dir.create(dir.name)
  file.name <- paste0(dir.name,"/",obj.name,".C.tsv")
  write_tsv(x, file = file.name, col_names=FALSE)
}


write.tsv.funct.D <- function(x){
  obj.name <- x$name %>% unique()
  dir.name <- snakemake@output[[1]]
  dir.create(dir.name)
  file.name <- paste0(dir.name,"/",obj.name,".D.tsv")
  write_tsv(x, file = file.name, col_names=FALSE)
}


lapply(bedpe.A, write.tsv.funct.A)

lapply(bedpe.B, write.tsv.funct.B)
lapply(bedpe.C, write.tsv.funct.C)
lapply(bedpe.D, write.tsv.funct.D)

# Saving the pre-processed data
# write_tsv(bed.file, snakemake@output[[1]], col_names=FALSE)

sink()
