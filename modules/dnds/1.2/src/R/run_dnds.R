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
    file = snakemake@input[["maf"]],
    sep = '\t',
    header = TRUE
)

maf = as.data.frame(maf)

if (grepl("38", maf$NCBI_Build[1])) {
    cat("Will run using the hg38 reference.\n")
    load(snakemake@params[[3]]) # This will load the RefCDS from the path at config
    use_this_reference <- RefCDS
}else{
    cat("Will run using the default reference.\n")
    use_this_reference <- "hg19"
}

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
        gene_list=target_genes$V1,
        max_muts_per_gene_per_sample=snakemake@params[[2]],
        refdb = use_this_reference
    )
}else{
    dndsout = dndscv(
        maf,
        max_muts_per_gene_per_sample=snakemake@params[[2]],
        refdb = use_this_reference
    )
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
