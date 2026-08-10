#!/usr/bin/env Rscript
# Parses HLA-HD's raw per-gene *.est.txt output into a clean {gene},{allele1},{allele2} CSV.
#
# This is module-owned code, NOT a copy or adaptation of MHC Hammer's own bin/hlahd_parse_output.R
# (which is off-limits to redistribute/modify under its Cancer Research Horizons academic-use
# licence -- see the licensing note near the top of mhc_hammer.smk). It was written for the HLA
# class II germline-typing path, which has no equivalent of upstream's IMGT-derived GTF (MHC
# Hammer's own Zenodo reference bundle, and therefore hlahd_parse_output.R's own tie-breaking
# logic, is class-I-only).
#
# Difference from hlahd_parse_output.R: when HLA-HD reports multiple ambiguous candidate alleles
# for a gene (comma-separated in its own est.txt output), hlahd_parse_output.R breaks ties by
# cross-referencing a GTF for which candidate has the most complete IMGT sequence representation.
# This script has no such GTF to reference, so it simply takes HLA-HD's own first-listed candidate
# per allele (HLA-HD already orders candidates by its own confidence) -- a deliberately simpler
# tie-break, not a like-for-like replacement.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument('--hlahd_folder', nargs = 1, required = TRUE,
                    help = 'Path to HLA-HD result directory containing {sample_id}_{gene}.est.txt files')
parser$add_argument('--sample_id', nargs = 1, required = TRUE)
parser$add_argument('--genes', nargs = '+', required = TRUE, help = 'Genes to parse, e.g. DRB1 DQA1 DQB1')
parser$add_argument('--output', nargs = 1, required = TRUE, help = 'Output CSV path')
args <- parser$parse_args()

not_typed <- function(gene) data.table(gene = gene, allele1 = "not typed", allele2 = "not typed")

result <- data.table(gene = character(), allele1 = character(), allele2 = character())

for (gene in args$genes) {
  est_file <- file.path(args$hlahd_folder, paste0(args$sample_id, "_", gene, ".est.txt"))

  if (!file.exists(est_file)) {
    result <- rbindlist(list(result, not_typed(gene)))
    next
  }

  lines <- readLines(est_file)

  if (length(lines) < 3 || lines[1] == "No candidate.") {
    result <- rbindlist(list(result, not_typed(gene)))
    next
  }

  fields <- strsplit(lines[3], "\t")[[1]]

  allele1_candidates <- strsplit(fields[1], ",")[[1]]
  allele1 <- trimws(allele1_candidates[1])

  if (length(fields) < 2 || fields[2] == "-") {
    # Homozygous (or only one allele group reported) -- both copies use the same call, matching
    # hlahd_parse_output.R's own handling of this case.
    allele2 <- allele1
  } else {
    allele2_candidates <- strsplit(fields[2], ",")[[1]]
    allele2 <- trimws(allele2_candidates[1])
  }

  result <- rbindlist(list(result, data.table(gene = gene, allele1 = allele1, allele2 = allele2)))
}

fwrite(result, args$output, col.names = FALSE)
cat("Saved to:", args$output, "\n")
