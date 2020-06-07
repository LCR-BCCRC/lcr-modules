#!/usr/bin/env Rscript

library(sequenza)
library(argparse)

# Parse arguments

parser <- ArgumentParser()

parser$add_argument("-i", "--input", 
                    help="input seqz file")
parser$add_argument("-o", "--outdir", 
                    help="output directory for sequenza output")
parser$add_argument("-n", "--name", 
                    help="sample name for file prefix")

args <- parser$parse_args()

seqz <- args$input
outdir <- args$outdir
samp <- args$name

# Run Sequenza

test <- sequenza.extract(seqz, verbose = FALSE)

CP <- sequenza.fit(test)

sequenza.results(sequenza.extract = test,
                 cp.table = CP, sample.id = samp,
                 out.dir=outdir)