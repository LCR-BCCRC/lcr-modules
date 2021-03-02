#!/usr/bin/env Rscript
#


suppressWarnings(suppressPackageStartupMessages({
  
  message("Loading packages...")
  
  required_packages <- c(
    "data.table", "dplyr", "tidyverse", "readr", "tidyr", "gtools"
  )
  
  installed_packages <- rownames(installed.packages())
  
  for (pkg in required_packages){
    if (!pkg %in% installed_packages) {
      message(paste0("Installing ", pkg, "..."))
      install.packages(pkg, repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
    }
    library(pkg, character.only = TRUE)
  }
}))


# Prepare data for CNVfilteR ----------------------------------------------
message("Loading data for cleanup after liftover...")

lifted.SEG <- fread(snakemake@input[[1]])
colnames(lifted.SEG) <- c("sample", "chrom", "start", "end", "LOH_flag", "log.ratio")

message("Preparing output...")
chromosomes <- c()
for (x in 1:22) {
  chromosomes[x] = paste0("chr",x)
  if (x==22) {
    chromosomes[23] = "chrX"
  }
}

output <- complete(lifted.SEG, tidyr::expand(lifted.SEG, crossing(sample), chrom = chromosomes))


output <- output %>% mutate(start=ifelse(is.na(start), paste0("100",rownames(.)), start)) %>%
  mutate(end=ifelse(is.na(end), paste0("100000",rownames(.)), end)) %>%
  mutate(LOH_flag=ifelse(is.na(LOH_flag), toString(0), LOH_flag)) %>%
  mutate(log.ratio=ifelse(is.na(log.ratio), 0.0, log.ratio))


# add new segment att at the end to allow for segment filling

new <- data.frame("sample" = output$sample[1],              
                  "chrom" = "chrX",
                  "start" = 156040865,
                  "end" = 156040875,
                  "LOH_flag" = 0,
                  "log.ratio" = 0.0)
output <- rbind(output, new)

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

output <- output[multi.mixedorder(output$chrom, output$start),]


# Output data ------------------------------------------------------
message("Writing final table to file")
write_tsv(output, snakemake@output[[1]])