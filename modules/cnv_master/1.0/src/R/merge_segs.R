#!/usr/bin/env Rscript
#


suppressWarnings(suppressPackageStartupMessages({

  message("Loading packages...")

  required_packages <- c(
    "readr", "dplyr", "tidyr", "purrr"
  )

  for (pkg in required_packages){
    library(pkg, character.only = TRUE)
  }
}))


# Read all individual seg files ----------------------------------------------
message("Loading data from individual seg files ...")

files <- snakemake@input[]

# read individual files and keep file paths in the column filename
data = lapply(files$seg_file, read_tsv, col_types = "ccddddd")


# strip file paths for the final seg file
output <- bind_rows(data) %>%
  distinct %>%
  as.data.frame

# this is the file path of all individual segs used in merging
contents <- data.frame(filename = files$seg_file)

# Output data ------------------------------------------------------
message("Writing final outputs ...")
write_tsv(output, snakemake@output[[1]])
write_tsv(contents, snakemake@output[[2]], col_names = FALSE)
