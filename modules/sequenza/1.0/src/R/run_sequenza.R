#!/usr/bin/env Rscript


# Load Packages -----------------------------------------------------------

remotes::install_github("ShixiangWang/copynumber", dependencies = FALSE, 
                        upgrade = "never")

library(sequenza)


# Parse Arguments ---------------------------------------------------------

args      <- commandArgs(trailingOnly=TRUE)
seqz_path <- args[1]
assembly  <- args[2]
chroms    <- args[3]
odir_path <- args[4]
num_cores <- args[5]


# Prepare variables -------------------------------------------------------

assembly <- readLines(assembly)

to_ucsc <- list(
  hg19 = "hg19",
  hg38 = "hg38",
  grch37 = "hg19",
  grch38 = "hg38"
)

stopifnot(assembly %in% names(to_ucsc))
assembly_ucsc <- to_ucsc[assembly]

chroms <- readLines(chroms)

if (is.null(num_cores)) num_cores <- 1
num_cores <- as.integer(num_cores)


# Load Data ---------------------------------------------------------------

seqz <- sequenza.extract(
  seqz_path,
  chromosome.list = chroms,
  assembly = assembly_ucsc
)


# Fit Cellularity and Ploidy ----------------------------------------------

chrx_ratio <- median(seqz$depths$norm$normal$X$mean, na.rm = TRUE)
is_female <- chrx_ratio >= 0.75
sex <- ifelse(is_female, "female", "male")
message(paste0("Inferring this patient as ", sex, "."))

fit <- sequenza.fit(
  seqz,
  female = is_female,
  mc.cores = num_cores
)


# Output Results ----------------------------------------------------------

res <- sequenza.results(
  seqz,
  fit,
  sample.id = "sequenza",
  out.dir = odir_path,
  female = is_female,
  chromosome.list = chroms
)
