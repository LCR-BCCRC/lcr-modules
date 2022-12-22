#!/usr/bin/env Rscript
#

log <- file(snakemake@log[[1]], open="wt")
sink(log)

suppressWarnings(
    suppressPackageStartupMessages(
        library(readr)
    )
)


# Read in the matrix
ge_matrix <- suppressMessages(
    read_tsv(
        snakemake@input[[1]]
    )
)
ge_matrix <- as.data.frame(
    ge_matrix
)

# Drop some columns if they are present in the matrix
if("gene_id" %in% colnames(ge_matrix)){
    ge_matrix <- ge_matrix[ , !names(ge_matrix) %in% "gene_id"]
}

if("ensembl_gene_id" %in% colnames(ge_matrix)){
    ge_matrix <- ge_matrix[ , !names(ge_matrix) %in% "ensembl_gene_id"]
}

# Read in the annotations
annotations <- suppressMessages(
    read_tsv(
        snakemake@input[[2]]
    )
)
annotations <- as.data.frame(
    annotations
)

# Create column ID if sample_id is present
if("sample_id" %in% colnames(annotations)){
    annotations$ID <- annotations$sample_id
}

# Make sure there are annotations for every sample
annotations <- subset(
    annotations,
    ID %in% colnames(ge_matrix)
)

# Subset matrix to these samples with annotations
ge_matrix <- ge_matrix[ , names(ge_matrix) %in% c(names(ge_matrix)[1], annotations$ID)]

# Drop duplicated gene names and those with NAs
ge_matrix <- ge_matrix[!duplicated(ge_matrix[,1]),]
ge_matrix <- ge_matrix[!is.na(ge_matrix[,1]),]

# Modify sample ids because ecotyper cannot handle those with
# special characters (. or -) and those starting from number

generate_names <- function(
    incoming_data
){
    # Drop all special characters
    processed_data <- gsub(
        '_|-|\\.',
        '',
        incoming_data
        )

    # Add string so it does not start from number
    processed_data <- paste0(
        "sample",
        processed_data
    )

    return(processed_data)
}

# Creating mapping data
mapping <- data.frame(
    sample_id = annotations$ID,
    mapping = generate_names(
        annotations$ID
    )
)

# Generate names properly handled by ecotyper
annotations$ID <- generate_names(
        annotations$ID
    )

colnames(ge_matrix) <- generate_names(
        colnames(ge_matrix)
    )

# Saving the pre-processed data and mapping
write_tsv(mapping, snakemake@output[[1]])
write_tsv(ge_matrix, snakemake@output[[2]])
write_tsv(annotations, snakemake@output[[3]])

sink()
