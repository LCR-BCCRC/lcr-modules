#!/usr/bin/env Rscript

# CNAqc CCF estimation script
# Inputs: Battenberg subclones (filled), cellularity_ploidy, MAF
# Output: per-mutation CCF TSV (all columns returned by CNAqc::CCF())

suppressPackageStartupMessages({
  library(optparse)
})

# ---- Parse arguments ----

option_list <- list(
  make_option("--subclones",         type = "character", help = "Battenberg subclones.txt file (filled)"),
  make_option("--cellularity_ploidy",type = "character", help = "Battenberg cellularity_ploidy.txt file"),
  make_option("--maf",               type = "character", help = "MAF file with somatic SNVs"),
  make_option("--tumour_id",         type = "character", help = "Tumour sample ID"),
  make_option("--out_ccf",           type = "character", help = "Output TSV file for CCF estimates"),
  make_option("--ref",               type = "character", default = "GRCh38",
              help = "CNAqc genome reference string: hg19, GRCh37, hg38, or GRCh38 [default: %default]"),
  make_option("--min_depth",         type = "integer",   default = 10L,
              help = "Minimum read depth for variants [default: %default]"),
  make_option("--min_muts_per_segment", type = "integer", default = 10L,
              help = "Minimum mutations per segment to attempt peak analysis [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

for (arg in c("subclones", "cellularity_ploidy", "maf", "tumour_id", "out_ccf")) {
  if (is.null(opt[[arg]])) stop(paste("Required argument missing:", arg))
}

suppressPackageStartupMessages({
  library(CNAqc)
  library(dplyr)
  library(data.table)
})

write_empty <- function(path, tumour_id) {
  write.table(
    data.frame(tumour_id = tumour_id),
    path, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# ---- Cellularity / ploidy ----

cp     <- read.table(opt$cellularity_ploidy, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
purity <- as.numeric(cp$cellularity[1])
ploidy <- as.numeric(cp$ploidy[1])
message(sprintf("Purity: %.4f  Ploidy: %.4f", purity, ploidy))

# ---- Battenberg subclones ----

sub         <- fread(opt$subclones, sep = "\t", data.table = FALSE)
colnames(sub) <- tolower(colnames(sub))

required_sub_cols <- c("chr", "startpos", "endpos", "nmaj1_a", "nmin1_a", "frac1_a")
missing_cols <- required_sub_cols[!required_sub_cols %in% colnames(sub)]
if (length(missing_cols) > 0) {
  stop(paste("Subclones file missing required columns:", paste(missing_cols, collapse = ", ")))
}

cna <- data.frame(
  chr     = sub$chr,
  from    = as.integer(sub$startpos),
  to      = as.integer(sub$endpos),
  Major   = as.integer(sub$nmaj1_a),
  minor   = as.integer(sub$nmin1_a),
  CCF     = as.numeric(sub$frac1_a),
  Major_2 = if ("nmaj2_a" %in% colnames(sub)) as.integer(sub$nmaj2_a) else NA_integer_,
  minor_2 = if ("nmin2_a" %in% colnames(sub)) as.integer(sub$nmin2_a) else NA_integer_,
  stringsAsFactors = FALSE
)

if (!all(grepl("^chr", cna$chr))) cna$chr <- paste0("chr", cna$chr)

cna <- cna[!is.na(cna$Major) & !is.na(cna$minor) & (cna$to - cna$from) > 1L, ]
message(sprintf("Loaded %d CNA segments.", nrow(cna)))

# ---- MAF / SNVs ----

skip_n  <- sum(grepl("^#", readLines(opt$maf, n = 200)))
maf_raw <- fread(opt$maf, sep = "\t", skip = skip_n, data.table = FALSE)
colnames(maf_raw) <- tolower(colnames(maf_raw))

find_col <- function(candidates) {
  m <- candidates[candidates %in% colnames(maf_raw)]
  if (length(m) == 0) return(NULL)
  m[1]
}

chr_col  <- find_col(c("chromosome", "chrom"))
from_col <- find_col(c("start_position", "start"))
ref_col  <- find_col(c("reference_allele"))
alt_col  <- find_col(c("tumor_seq_allele2", "alt_allele"))
dp_col   <- find_col(c("t_depth"))
nv_col   <- find_col(c("t_alt_count"))

for (col_name in c("chr_col", "from_col", "ref_col", "alt_col", "dp_col", "nv_col")) {
  if (is.null(get(col_name))) stop(paste("Cannot find required MAF column for:", col_name))
}

snvs <- data.frame(
  chr  = maf_raw[[chr_col]],
  from = as.integer(maf_raw[[from_col]]),
  to   = as.integer(maf_raw[[from_col]]) + 1L,
  ref  = maf_raw[[ref_col]],
  alt  = maf_raw[[alt_col]],
  DP   = as.integer(maf_raw[[dp_col]]),
  NV   = as.integer(maf_raw[[nv_col]]),
  stringsAsFactors = FALSE
)

snvs$VAF <- snvs$NV / snvs$DP
snvs <- snvs[
  !is.na(snvs$DP) & snvs$DP >= opt$min_depth &
    nchar(snvs$ref) == 1L & nchar(snvs$alt) == 1L, ]

if (!all(grepl("^chr", snvs$chr))) snvs$chr <- paste0("chr", snvs$chr)
message(sprintf("Loaded %d SNVs (depth >= %d).", nrow(snvs), opt$min_depth))

# ---- CNAqc init ----

x <- CNAqc::init(mutations = snvs, cna = cna, purity = purity, ref = opt$ref)

# ---- Peak analysis (required before compute_CCF) ----

if (nrow(snvs) < opt$min_muts_per_segment) {
  message("Too few mutations for peak analysis — writing empty CCF file.")
  write_empty(opt$out_ccf, opt$tumour_id)
  quit(save = "no", status = 0)
}

peaks_ok <- tryCatch({
  x <- CNAqc::analyze_peaks(x, matching_strategy = "closest")
  TRUE
}, error = function(e) {
  message(paste("Peak analysis failed:", conditionMessage(e)))
  FALSE
})

if (!peaks_ok) {
  message("Writing empty CCF file.")
  write_empty(opt$out_ccf, opt$tumour_id)
  quit(save = "no", status = 0)
}

# ---- CCF estimation ----

ccf_ok <- tryCatch({
  x <- CNAqc::compute_CCF(x)
  TRUE
}, error = function(e) {
  message(paste("compute_CCF failed:", conditionMessage(e)))
  FALSE
})

if (!ccf_ok) {
  message("Writing empty CCF file.")
  write_empty(opt$out_ccf, opt$tumour_id)
  quit(save = "no", status = 0)
}

# ---- Export ----

ccf_df <- tryCatch(
  CNAqc::CCF(x),
  error = function(e) {
    message(paste("CCF() getter failed:", conditionMessage(e)))
    NULL
  }
)

if (is.null(ccf_df) || nrow(ccf_df) == 0) {
  message("No CCF estimates available — writing empty CCF file.")
  write_empty(opt$out_ccf, opt$tumour_id)
} else {
  ccf_df <- cbind(tumour_id = opt$tumour_id, ccf_df)
  write.table(ccf_df, opt$out_ccf, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("Wrote %d CCF estimates to %s", nrow(ccf_df), opt$out_ccf))
}

message("Done.")
