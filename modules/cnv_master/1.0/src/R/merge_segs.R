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
data =
tibble(filename = files$seg_file) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(files$seg_file, col_types = cols())) # a new data column
  ) %>%
  unnest(cols = c(file_contents))

# strip file paths for the final seg file
output <- data %>%
  select(-filename) %>%
  distinct %>%
  as.data.frame

# this is the file path of all individual segs used in merging
contents <- data %>%
  select(filename) %>%
  unique %>%
  as.data.frame

# Output data ------------------------------------------------------
message("Writing final outputs ...")
write_tsv(output, snakemake@output[[1]])
write_tsv(contents, snakemake@output[[2]], col_names = FALSE)
