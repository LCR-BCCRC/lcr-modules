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
    library(fishHook)
    library(gTrack)
    library(rtracklayer)
  })
)
#load maf and convert to GRange
maf = dt2gr(fread(snakemake@input[[2]]))

#split maf to numerous tiles
tiles = gr.tile(seqinfo(maf), snakemake@params[[1]])

events = maf %Q% (Variant_Classification != 'Silent')

if(file.exists(snakemake@params[[2]])){
  message("Running FishHook with Coveriate...")
  fish_tiles = Fish(hypotheses = tiles, 
                    events = events, 
                    eligible = eligible, 
                    idcol = 'Tumor_Sample_Barcode',
                    use_local_mut_density=TRUE)
  
  coveriate_data = gr.sub(import(snakemake@params[[2]]), 'chr', "") ## import from bed then gUtils::gr.sub to strip 'chr' identifier
  coveriate = Cov(coveriate_data, name = 'coveriate')
  fish_tiles$covariates = c(coveriate)
  
}else{
  message("Running FishHook without Coveriate...")
  fish_tiles = Fish(hypotheses = tiles, 
                    events = events, 
                    eligible = eligible, 
                    idcol = 'Tumor_Sample_Barcode',
                    use_local_mut_density=TRUE)
  
}

fish_tiles$score()

result = fish_tiles$res %Q% which(p<=0.05)

df.result = gr2dt(result)

write.table(
  df.result,
  file=snakemake@output[[1]],
  quote=FALSE,
  sep='\t',
  row.names = FALSE)



