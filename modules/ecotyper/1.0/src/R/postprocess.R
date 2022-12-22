#!/usr/bin/env Rscript
#

suppressWarnings(
    suppressPackageStartupMessages(
        library(readr)
    )
)

# Get the mapping to sample ids
mapping <- read_tsv(
    normalizePath(
        snakemake@input[[1]]
    )
)

# Sample file to get the path to the output directory
file <- snakemake@input[[2]]

# Get the second level direstory of all outputs
dir <- gsub(
    "mapped_ge_matrix.*",
    "mapped_ge_matrix",
    file
)

# Helper function to map all column names
map_names <- function(file_path){
    incoming_file <- read.delim(file_path)

    arranged_mapping <- mapping[match(colnames(incoming_file), mapping$mapping),]

    colnames(incoming_file) <- arranged_mapping$sample_id

    write.table(incoming_file, file_path, sep = "\t")
}

# Helper function to map all row names and IDs
map_ids <- function(file_path){
    incoming_file <- read.delim(file_path)

    arranged_mapping <- mapping[match(incoming_file$ID, mapping$mapping),]

    rownames(incoming_file) <- arranged_mapping$sample_id

    incoming_file$ID <- arranged_mapping$sample_id

    write.table(incoming_file, file_path, sep = "\t")
}

# Converting all wide files
files <- list.files(dir,
                    recursive = TRUE,
                    pattern = "heatmap_data|abundance",
                    full.names = TRUE)

lapply(
    files,
    map_names
)

# Converting all long files
files <- list.files(dir,
                    recursive = TRUE,
                    pattern = "heatmap_top_ann|assignment.txt",
                    full.names = TRUE)

lapply(
    files,
    map_ids
)

# Writing sentinnel file to indicate the end of snakemake rule
cat(
    NULL,
    file = snakemake@output[[1]]
)