#! /usr/bin/env Rscript

script_path <- normalizePath(dirname(unlist(strsplit(commandArgs()[grepl('--file', commandArgs())], split = '='))[2]))

args <- commandArgs(trailingOnly = TRUE)

threads  <- args[1]
nrml.dir <- args[2]
tbl.out  <- args[3]
pon.path <- args[4]
par.path <- args[5]
build    <- args[6]

library(data.table)
library(dryclean)
library(parallel)

nrml.dir.abs <- normalizePath(nrml.dir)
nrml.rds <- Sys.glob(file.path(nrml.dir.abs, "*.rds"))
nrml.smp <- gsub('.cov.rds', '', basename(nrml.rds))

nrml.tbl <- data.table(
    sample = nrml.smp,
    normal_cov = nrml.rds
)

saveRDS(nrml.tbl, tbl.out)

source(file.path(script_path, 'prepare_detergent.R'))

detergent = prepare_detergent(
    normal.table.path = tbl.out,
    path.to.save = dirname(pon.path),
    num.cores = threads,
    build = build,
    PAR.file = par.path,
    save.pon = TRUE,
    choose.by.clustering = TRUE
    #use.all = TRUE
)


