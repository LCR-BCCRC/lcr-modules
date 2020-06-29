#!/usr/bin/env Rscript


# Load Packages -----------------------------------------------------------

remotes::install_github("ShixiangWang/copynumber", dependencies = FALSE)

library(sequenza)


# Parse Arguments ---------------------------------------------------------

args      <- commandArgs(trailingOnly=TRUE)
seqz_path <- args[[1]]
assembly  <- args[[2]]
chroms    <- args[[3]]
x_chrom   <- args[[4]]
odir_path <- args[[5]]
num_cores <- args[[6]]


# Prepare variables -------------------------------------------------------

assembly <- readLines(assembly)

to_ucsc <- list(
  hg19 = "hg19",
  hg38 = "hg38",
  grch37 = "hg19",
  grch38 = "hg38"
)

stopifnot(assembly %in% names(to_ucsc))
assembly_ucsc <- to_ucsc[[assembly]]
assembly_ucsc <- unname(assembly_ucsc)

chroms <- readLines(chroms)

x_chrom <- readLines(x_chrom)

if (is.null(num_cores)) num_cores <- 1
num_cores <- as.integer(num_cores)


# Load Data ---------------------------------------------------------------

seqz <- sequenza.extract(
  seqz_path,
  chromosome.list = chroms,
  assembly = assembly_ucsc
)


# Fit Cellularity and Ploidy ----------------------------------------------

num_positions <- sum(seqz$depths$norm$normal[[x_chrom]]$N, na.rm = TRUE)
if (num_positions < 1e6) {
  warn <- paste0("There are fewer than 1 million positions with depth data on the ",
                 "X chromosome. This usually means the input data was based on ",
                 "exomes, which will fail the sex inferrence step.")
  warning(warn)
}

chrx_ratio <- median(seqz$depths$norm$normal[[x_chrom]]$mean, na.rm = TRUE)
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
