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

my_merge_function <- function(path) {
    col_count <- count.fields(path, sep = "\t")[1]
    if(col_count == 9){
        # ID      chrom   start   end     module  LOH_flag        CN      log.ratio       dummy_segment
        incoming_data <- suppressMessages(
        suppressWarnings(
            read_tsv(
                path,
                col_types = "cciiciidi",
                progress = FALSE
                )
            )
        )
    }else if(col_count == 7){
    # ID chrom start end num.mark seg.mean dummy_segment
    # ID chrom loc.start loc.end num.mark seg.mean dummy_segment
    # ID chrom start end  LOH_flag  log.ratio  dummy_segment
        incoming_data <- suppressMessages(
            suppressWarnings(
                read_tsv(
                    path,
                    col_types = "cciiidi",
                    progress = FALSE
                )
            )
        )

        colnames(incoming_data) <- gsub(
            "loc.",
            "",
            colnames(incoming_data)
        )

        if ("seg.mean" %in% colnames(incoming_data)) {
        # ID chrom start end num.mark seg.mean dummy_segment
            incoming_data <- rename(
                incoming_data,
                log.ratio = seg.mean
            ) %>%
            select(-num.mark) %>%
            mutate(
                LOH_flag = NA,
                module = NA,
                CN = NA
            )
        }else{
        # ID chrom start end  LOH_flag  log.ratio  dummy_segment
            incoming_data <- incoming_data %>%
                mutate(
                    module = NA,
                    CN = NA
                )
        }
        incoming_data <- incoming_data %>%
            select(ID,chrom,start,end,module,LOH_flag,CN,log.ratio,dummy_segment)
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
