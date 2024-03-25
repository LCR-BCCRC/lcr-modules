### Run R package fishHook to indentify statistical enrichment or depletion of somatic mutations in arbitrary (sets of) genomic intervals ###

#!/usr/bin/env Rscript

# Notes:
#   This script is intended for use with the snakefile for fishHook in lcr-modules.
#   It uses the snakemake objects for the input maf file and paramters.
#   Additionaly, it can support the genes of targeted panel to restrict SMG search.

log <- file(snakemake@log[[1]], open="wt")
sink(log ,type = "output")
sink(log, type = "message")

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressWarnings(
  suppressPackageStartupMessages({
    library(fishHook)
    library(gTrack)
    library(rtracklayer)
  })
)

# Read in maf and convert to GRange
maf = gr.sub(dt2gr(fread(snakemake@input[[2]])))

if(!snakemake@params[[3]]){
  message("Excluding Silent Mutations ...")
  events = maf %Q% (Variant_Classification != 'Silent')
}else{
  events = maf
}

# Use tile mode or gene list mode
if(is.null(snakemake@params[[1]])){
  message("Running FishHook with Gene List...")
  genes = gr.sub(import(snakemake@params[[4]]))

  if(snakemake@params[[5]]){
    message("Subsetting Gene List for Protein Coding Gene Only ...")
    genes = genes %Q% (gene_type == 'protein_coding')
  }
  fish = Fish(hypotheses = genes,
                    events = events,
                    idcol = 'Tumor_Sample_Barcode',
                    use_local_mut_density=TRUE)
}else{
  message("Running FishHook with Tiles ...")
  # Split maf to tiles
  tiles = gr.tile(seqinfo(maf), snakemake@params[[1]])

  fish = Fish(hypotheses = tiles,
                    events = events,
                    idcol = 'Tumor_Sample_Barcode',
                    use_local_mut_density=TRUE)
}

# If user provide covariates file
if(file.exists(snakemake@params[[2]])){
  message("Running FishHook with Coveriate...")
  covariates_data = gr.sub(import(snakemake@params[[2]]), 'chr', "") ## import from bed then gUtils::gr.sub to strip 'chr' identifier
  covariate = Cov(covariates_data, name = 'covariate')
  fish$covariates = c(covariate)

}else{
  message("Running FishHook without Covariates...")
}

fish$score()

result = fish$res %Q% which(p<=0.05) %Q% order(p) # subset result to include only region with p<0.05

df.result = gr2dt(result)

write.table(
  df.result,
  file=snakemake@output[[1]],
  quote=FALSE,
  sep='\t',
  row.names = FALSE)

sink()

