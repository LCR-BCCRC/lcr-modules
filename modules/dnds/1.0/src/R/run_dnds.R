### Run R package dNdS to determine the significantly mutated genes ###
# Currently supports only hg19-based maf files
# but future versions can be adopted to work with other genome builds.

#!/usr/bin/env Rscript
#
# Notes:
#   This script is intended for use with the snakefile for dNdS in LCR-modules.
#   It uses the snakemake objects to find maf file.
#   Additionaly, it can support the genes of targeted panel to restrict SMG search.
#

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressWarnings(
suppressPackageStartupMessages({
    library(dndscv)
})
)

# Read the maf file
maf = read.table(
    file = snakemake@input[[2]],
    sep = '\t',
    header = TRUE
)

maf = as.data.frame(maf)

# Run the dNdS
dndsout = dndscv(maf)

# Sig genes
sig_genes = as.data.frame(dndsout$sel_cv)

write_tsv(sig_genes,
          snakemake@output[[1]]
)

# Annotated mutations
annotmuts = as.data.frame(dndsout$annotmuts)

write_tsv(sig_genes,
          snakemake@output[[2]]
)
