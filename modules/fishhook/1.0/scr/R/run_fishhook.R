### Run R package dNdS to determine the significantly mutated genes ###

#!/usr/bin/env Rscript
#
# Notes:
#   This script is intended for use with the snakefile for fishHook in LCR-modules.
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