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

# This function will handle discrepancy in the output from
# CNVkit and Pure CN and will harmonize the colnames, output format

my_merge_function <- function(path) {
    incoming_data <- suppressMessages(
        suppressWarnings(
            read_tsv(
                path,
                col_types = "ccddddd"
            )
        )
    )

    colnames(incoming_data) <- gsub(
        "loc.",
        "",
        colnames(incoming_data)
    )

    if ("seg.mean" %in% colnames(incoming_data)) {
        incoming_data <- rename(
            incoming_data,
            log.ratio = seg.mean
        )
    }

    if ("num.mark" %in% colnames(incoming_data)) {
        incoming_data <- select(
            incoming_data,
            -num.mark
        )

        incoming_data <- mutate(
            incoming_data,
            LOH_flag = NA,
            .before = "log.ratio"
        )
    }

    return(incoming_data)
}

data <- lapply(
  files$seg_file,
  my_merge_function
)

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
