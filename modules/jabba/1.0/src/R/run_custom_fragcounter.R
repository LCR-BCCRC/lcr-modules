#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
    library(GenomicRanges)
    library(gUtils)
    library(rtracklayer)
    library(GenomeInfoDb)
    library(data.table)
    library(stats)
    library(skidb)
    library(Biostrings)
    library(parallel)
    library(grDevices)
    library(IRanges)
    library(DNAcopy)
    library(graphics)
    library(utils)
    library(bamUtils)
    library(Rsamtools)
}))

script_path <- normalizePath(dirname(unlist(strsplit(commandArgs()[grepl('--file', commandArgs())], split = '='))[2]))

args <- commandArgs(trailingOnly = TRUE)

bam    <- args[1]
window <- args[2]
gcmap  <- args[3]
outdir <- args[4]

source(file.path(script_path, 'custom_fragcounter.R'))

fragCounter(bam = bam, 
            window = as.numeric(window),
            gc.rds.dir = gcmap,
            map.rds.dir = gcmap,
            outdir = outdir)
