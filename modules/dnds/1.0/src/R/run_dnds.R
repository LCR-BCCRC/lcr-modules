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
# should we restrict to target genes?
if(file.exists(snakemake@params[[1]])){
    print(paste("Restricting analysis to target genes specified in", snakemake@params[[1]]))
    # Read the target genes file
    target_genes = read.table(
        file = snakemake@params[[1]],
        sep = '\t',
        header = FALSE
    )
    dndsout = dndscv(
        maf,
        gene_list=target_genes$V1
    )
}else{
    dndsout = dndscv(maf)
}

# Sig genes
sig_genes = as.data.frame(dndsout$sel_cv)

write.table(
    sig_genes,
    file=snakemake@output[[1]],
    quote=FALSE,
    sep='\t',
    row.names = FALSE)

# Annotated mutations
annotmuts = as.data.frame(dndsout$annotmuts)

write.table(
    annotmuts,
    file=snakemake@output[[2]],
    quote=FALSE,
    sep='\t',
    row.names = FALSE)
