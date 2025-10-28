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

save.image("/projects/rmorin/projects/gambl-repos/gambl-lhilton/results/all_the_things/dlbclass-1.0/01-prepare_inputs/gambl-all--2025-10/seg2gsm.RData")

# Load sample set file
samples <- read.delim(samples, sep = "\t")
samples <- samples[[1]]

# Load SV file
svs <- lapply(sv_files, read.delim, sep = "\t")
svs <- do.call(rbind, svs)
svs <- unique(svs)

# Confirm correct dimensions
if (length(colnames(svs)) != 4) {
    stop("SV file must have 4 columns with sample_id and MYC, BCL2, and BCL6 SV status. ")
}

# Confirm sample_id column
if (!"sample_id" %in% colnames(svs)) {
    stop("SV file must have a sample_id column. ")
}

# Confirm MYC, BCL2, BCL6 columns
if (!all(sapply(c("MYC", "BCL2", "BCL6"), function(x) any(grepl(x, colnames(svs), ignore.case = TRUE))))) {
    stop("SV file must have columns containing 'MYC', 'BCL2', and 'BCL6' as substrings.")
}

# Convert column names to match gsm requirements
colnames(svs) <- gsub(".*(BCL2|BCL6|MYC).*", "\\1", toupper(colnames(svs)))
colnames(svs) <- paste0("SV.", colnames(svs))
colnames(svs) <- gsub("SV.SAMPLE_ID", "sample_id", colnames(svs))

# Convert SV values to GSM values (0 for negative, 3 for positive)
svs$SV.BCL6 <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.BCL6))
svs$SV.BCL2 <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.BCL2))
svs$SV.MYC <- gsub(neg_value, 0, gsub(pos_value, 3, svs$SV.MYC))

# Pivot to feature rows and sample columns
rownames(svs) <- svs$sample_id
svs <- svs[, c("SV.BCL6", "SV.BCL2", "SV.MYC")]
samples <- samples[samples %in% rownames(svs)]
svs <- data.frame(t(svs[samples, ]), check.names = FALSE)
svs$classifier_name <- rownames(svs)
svs <- svs[, c("classifier_name", samples)]

# Write to file
write.table(svs, snakemake@output[["sv_gsm"]], row.names = FALSE, quote = FALSE, sep = "\t")

sink()
