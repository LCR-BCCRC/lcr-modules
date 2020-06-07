#!/usr/bin/env Rscript


# Load Packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(sequenza)
})


# Parse Arguments ---------------------------------------------------------

args      <- commandArgs(trailingOnly=TRUE)
seqz_path <- args[1] 
odir_path <- args[2]
sex       <- args[3] 
num_cores <- args[4] 
sample_id <- args[5]
assembly  <- args[6] 

num_cores <- as.integer(num_cores)

print(args)
# Prepare variables -------------------------------------------------------

is_female  <- sex %in% c("female", "Female", "FEMALE", "f", "F")

#there must be a more elegant way to handle this that is lcr-modules friendly
if(assembly == "hg19" || assembly == "grch37"){
  chroms     <- c(seq(1,22),"X")
  if(assembly == "grch37"){
	assembly = "hg19";
  }
} else{
  chroms     <- paste0("chr",c(seq(1,22), "X"))
}

# Load Data ---------------------------------------------------------------

seqz <- sequenza.extract(
  seqz_path,
  chromosome.list = chroms,
  assembly = assembly
)


# Fit Cellularity and Plloidy ---------------------------------------------

fit <- sequenza.fit(
  seqz,
  female = is_female,
  mc.cores = num_cores
)


# Output Results ----------------------------------------------------------

res <- sequenza.results(
  seqz,
  fit,
  sample.id = sample_id,
  out.dir = odir_path,
  female = is_female,
  chromosome.list = chroms
)
