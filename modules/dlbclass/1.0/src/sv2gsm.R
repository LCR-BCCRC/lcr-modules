#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load Snakemake values
samples <- snakemake@input[[2]]
pos_value <- snakemake@params[["pos_value"]]
neg_value <- snakemake@params[["neg_value"]]
sv_files <- snakemake@input[["sv"]]

# samples <- "/home/lhilton/repos/DLBCLass/DLBclass-tool/DLBCL_GenomeCanada/DLBCL_GenomeCanada_sampleset.tsv"
# pos_value <- "POS"
# neg_value <- "NEG"

# Load sample set file
samples <- read.delim(samples, sep = "\t")
samples <- samples[[1]]

# Load SV file
svs <- lapply(sv_files, read.delim, sep = "\t")
svs <- do.call(rbind, svs)
svs <- svs[!duplicated(svs[, 1]), ]

# Confirm correct dimensions
if (length(colnames(svs)) != 4) {
    stop("SV file must have 4 columns with sample_id and MYC, BCL2, and BCL6 SV status. ")
}

# Convert column names to match gsm requirements
colnames(svs) <- gsub(".*(BCL2|BCL6|MYC).*", "\\1", toupper(colnames(svs)))
colnames(svs) <- paste0("SV.", colnames(svs))
colnames(svs)[1] <- "sample_id"

# Convert SV values to GSM values (0 for negative, 3 for positive)
svs$SV.BCL6 <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.BCL6))
svs$SV.BCL2 <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.BCL2))
svs$SV.MYC <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.MYC))

# Match df to input samples
rownames(svs) <- svs$sample_id
svs <- svs[, 2:4]
samples_keep <- samples[samples %in% rownames(svs)]
message(paste0("Number of samples in SV file: ", length(samples_keep)))

samples_missing <- samples[!samples %in% samples_keep]
message(paste0("Number of samples missing from SV file: ", length(samples_missing)))

# Make a matrix for missing samples
sv_empty <- data.frame(matrix(0, ncol = 1, nrow = 3))
rownames(sv_empty) <- colnames(svs)
if (length(samples_missing) > 0) {
    sv_empty <- data.frame(matrix(0, ncol = length(samples_missing), nrow = 3))
    colnames(sv_empty) <- samples_missing
    rownames(sv_empty) <- colnames(svs)
}
message(paste0("Dimensions of sv_empty: ", paste0(dim(sv_empty), collapse = "x")))

# Convert non-empty SV DF to GSM format
if (length(samples_keep) > 0) {
    svs <- data.frame(t(svs[samples_keep, ]), check.names = FALSE)
    svs <- cbind(sv_empty, svs)
} else {
    (svs <- sv_empty)
}
message(paste0("Dimensions of svs: ", paste0(dim(svs), collapse = "x")))

# Final formatting steps
svs$classifier_name <- rownames(svs)
svs <- svs[, c("classifier_name", samples)]


# Write to file
write.table(svs, snakemake@output[["sv_gsm"]], row.names = FALSE, quote = FALSE, sep = "\t")

sink()
