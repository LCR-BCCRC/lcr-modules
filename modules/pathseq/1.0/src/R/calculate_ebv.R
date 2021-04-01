#!/usr/bin/env Rscript
#

# This script will use number of mapped reads from flagstat and Pathseq output to calculate the EBV status of the sample
# based on the cut-off values that are provided by user in config file. It assumes 3 types of possible EBV status, being
# EBV-positive, EBV-intermediate, and EBV-negative.


suppressWarnings(suppressPackageStartupMessages({
  
  message("Loading packages...")
  
  required_packages <- c(
    "data.table", "tidyverse"
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


# this is used to assign EBV-status below
CUTOFF <-(str_split_fixed(sub("\\ +,*", "",snakemake@params[[1]]), ",", 2))
CUTOFF <- c(-Inf, as.numeric(CUTOFF[1]), as.numeric(CUTOFF[2]),  Inf)
# I am also saving the cutoff values separated by pipe character to provide this information in the output file
cutoff <- paste0(format(CUTOFF[2], scientific = FALSE), "|", format(CUTOFF[3], scientific = FALSE))


# read flagstat file and get the number of mapped reads
FLAGSTAT <- readLines(snakemake@input[[2]])
mapped_reads <- as.numeric(sub("\\ +.*", "",FLAGSTAT[5]))


# Load the Pathseq report
DATA <- read.delim(file=snakemake@input[[1]], sep = '\t', header = TRUE)


# read the pathseq scores files and get the number of EBV latency I reads
# for some files, there are no reads and EBV is missing

# first, initialize variables
PATHSEQ.RESULTS <- data.frame(matrix(ncol=10,nrow=0))
colnames(PATHSEQ.RESULTS) <- colnames(DATA)

# Search the coordinates of the cell with EBV report
COORDINATES <- which (DATA=="Human_gammaherpesvirus_4", arr.ind = TRUE)
# If no EBV in report file, add empty row to the DF with Pathseq results and replace the NA value with 0
if (all(is.na(COORDINATES))){
  PATHSEQ.RESULTS[1,]<-0
}else{
  # If EBV is found in the report file, append the whole row to the DF with Pathseq results    
  PATHSEQ.RESULTS <- rbind(PATHSEQ.RESULTS, DATA[COORDINATES[1,1],])
}

# Assemble the output table
PATHSEQ.RESULTS <- cbind(PATHSEQ.RESULTS, snakemake@wildcards$sample_id, mapped_reads) %>%
  `colnames<-`(c(colnames(PATHSEQ.RESULTS), "sample_id", "mapped_reads")) %>%
  # rename column with reads and calsulate fraction of ebv reads
  mutate(ebv_reads=unambiguous,
         ebv_reads_fraction=ebv_reads/mapped_reads) %>%
  # only select column that needed in the output
  select(sample_id, mapped_reads, ebv_reads, ebv_reads_fraction) %>%
  # save the information about cutoff values in final output
  cbind(.,cutoff)

# assign EBV status
PATHSEQ.RESULTS$ebv_status <- cut(PATHSEQ.RESULTS$ebv_reads_fraction, breaks=CUTOFF, label=c("EBV-Negative", "EBV-Intermediate", "EBV-Positive"))

# write final file
write_tsv(PATHSEQ.RESULTS, snakemake@output[[1]])
