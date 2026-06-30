#!/usr/bin/env Rscript

# CNAqc analysis script
# Inputs:
#   - Battenberg subclones (filled), cellularity_ploidy, MAF file
# Outputs:
#   - QC plots (PDF), metrics/score (TSV)

suppressPackageStartupMessages({
  library(optparse)
})

# ---- Parse arguments ----

option_list <- list(
  make_option(
    c("--subclones"),
    type = "character",
    help = "Battenberg subclones.txt file (filled)"
  ),
  make_option(
    c("--cellularity_ploidy"),
    type = "character",
    help = "Battenberg cellularity_ploidy.txt file"
  ),
  make_option(
    c("--maf"),
    type = "character",
    help = "MAF file with somatic SNVs"
  ),
  make_option(
    c("--tumour_id"),
    type = "character",
    help = "Tumour sample ID"
  ),
  make_option(
    c("--out_plot"),
    type = "character",
    help = "Output PDF file for QC plots"
  ),
  make_option(
    c("--out_metrics"),
    type = "character",
    help = "Output TSV file for QC metrics"
  ),
  make_option(
    c("--ref"),
    type = "character",
    default = "GRCh38",
    help = paste0(
      "CNAqc genome reference string: ",
      "hg19, GRCh37, hg38, or GRCh38 [default: %default]"
    )
  ),
  make_option(
    c("--min_depth"),
    type = "integer",
    default = 10L,
    help = "Minimum read depth for variants [default: %default]"
  ),
  make_option(
    c("--min_muts_per_segment"),
    type = "integer",
    default = 10L,
    help = "Minimum mutations per segment [default: %default]"
  ),
  make_option(
    c("--run_peak_analysis"),
    type = "logical",
    default = TRUE,
    help = "Run peak analysis [default: %default]"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = paste0(
    "Run CNAqc QC analysis using ",
    "Battenberg CNA and MAF SNVs"
  )
)
opt <- parse_args(opt_parser)

required_args <- c(
  "subclones", "cellularity_ploidy", "maf",
  "tumour_id", "out_plot", "out_metrics"
)
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    stop(paste("Required argument missing:", arg))
  }
}

valid_refs <- c("hg19", "GRCh37", "hg38", "GRCh38")
if (!opt$ref %in% valid_refs) {
  stop(paste(
    "--ref must be one of:", paste(valid_refs, collapse = ", ")
  ))
}

suppressPackageStartupMessages({
  library(CNAqc)
  library(dplyr)
  library(ggplot2)
  library(data.table)
})

message("Loading input files...")

# ---- Cellularity / ploidy ----

cp <- read.table(
  opt$cellularity_ploidy,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
purity <- as.numeric(cp$cellularity[1])
ploidy  <- as.numeric(cp$ploidy[1])
message(sprintf("Purity: %.4f  Ploidy: %.4f", purity, ploidy))

# ---- Battenberg subclones ----

sub <- fread(opt$subclones, sep = "\t", data.table = FALSE)
colnames(sub) <- tolower(colnames(sub))

required_sub_cols <- c(
  "chr", "startpos", "endpos",
  "nmaj1_a", "nmin1_a", "frac1_a"
)
missing_cols <- required_sub_cols[
  !required_sub_cols %in% colnames(sub)
]
if (length(missing_cols) > 0) {
  stop(paste(
    "Subclones file missing required columns:",
    paste(missing_cols, collapse = ", ")
  ))
}

# CNAqc expects: chr, from, to, Major, minor, CCF, Major_2, minor_2
# Major_2/minor_2 come from the second Battenberg clone (nmaj2_a/nmin2_a);
# filled with NA for purely clonal segments or absent columns.
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

# CNAqc requires chr-prefixed chromosome names
if (!all(grepl("^chr", cna$chr))) {
  cna$chr <- paste0("chr", cna$chr)
}

# Drop NA copy numbers or zero-length segments
cna <- cna[
  !is.na(cna$Major) &
    !is.na(cna$minor) &
    (cna$to - cna$from) > 1L,
]

message(sprintf(
  "Loaded %d CNA segments from subclones file.", nrow(cna)
))

# ---- MAF / SNVs ----

skip_n <- sum(grepl("^#", readLines(opt$maf, n = 200)))
maf_raw <- fread(
  opt$maf, sep = "\t", skip = skip_n, data.table = FALSE
)
colnames(maf_raw) <- tolower(colnames(maf_raw))
maf_cols <- colnames(maf_raw)

find_col <- function(candidates) {
  m <- candidates[candidates %in% maf_cols]
  if (length(m) == 0) return(NULL)
  m[1]
}

chr_col  <- find_col(c("chromosome", "chrom"))
from_col <- find_col(c("start_position", "start"))
ref_col  <- find_col(c("reference_allele"))
alt_col  <- find_col(c("tumor_seq_allele2", "alt_allele"))
dp_col   <- find_col(c("t_depth"))
nv_col   <- find_col(c("t_alt_count"))

for (col_name in c(
  "chr_col", "from_col", "ref_col",
  "alt_col", "dp_col", "nv_col"
)) {
  if (is.null(get(col_name))) {
    stop(paste(
      "Cannot find required MAF column for:", col_name
    ))
  }
}

snvs <- data.frame(
  chr = maf_raw[[chr_col]],
  from = as.integer(maf_raw[[from_col]]),
  to   = as.integer(maf_raw[[from_col]]) + 1L,
  ref  = maf_raw[[ref_col]],
  alt  = maf_raw[[alt_col]],
  DP   = as.integer(maf_raw[[dp_col]]),
  NV   = as.integer(maf_raw[[nv_col]]),
  stringsAsFactors = FALSE
)

snvs$VAF <- snvs$NV / snvs$DP

# Keep SNVs only (single-base substitutions) above depth threshold
snvs <- snvs[
  !is.na(snvs$DP) &
    snvs$DP >= opt$min_depth &
    nchar(snvs$ref) == 1L &
    nchar(snvs$alt) == 1L,
]

if (!all(grepl("^chr", snvs$chr))) {
  snvs$chr <- paste0("chr", snvs$chr)
}

message(sprintf(
  "Loaded %d SNVs (depth >= %d).", nrow(snvs), opt$min_depth
))

# ---- CNAqc init ----

message(sprintf(
  "Initialising CNAqc object (ref = %s)...", opt$ref
))
x <- CNAqc::init(
  mutations = snvs,
  cna       = cna,
  purity    = purity,
  ref       = opt$ref
)

# ---- Peak analysis ----

if (opt$run_peak_analysis &&
    nrow(snvs) >= opt$min_muts_per_segment) {
  message("Running peak analysis...")
  tryCatch(
    x <- CNAqc::analyze_peaks(
      x, matching_strategy = "closest"
    ),
    error = function(e) {
      message(paste("Peak analysis failed:", conditionMessage(e)))
    }
  )
}

# ---- CNAqc QC metrics ----

message("Computing CNAqc QC metrics...")
pass_pct <- tryCatch(
  CNAqc::get_PASS_percentage(x),
  error = function(e) {
    message(paste("get_PASS_percentage failed:", conditionMessage(e)))
    NULL
  }
)

score_df <- data.frame(
  tumour_id      = opt$tumour_id,
  ref            = opt$ref,
  purity         = purity,
  ploidy         = ploidy,
  n_snvs         = nrow(snvs),
  n_cna_segments = nrow(cna),
  pass_pct       = if (!is.null(pass_pct) && nrow(pass_pct) > 0)
                     pass_pct$percentage[pass_pct$QC_PASS %in% c("PASS", TRUE)][1]
                   else NA_real_,
  stringsAsFactors = FALSE
)

# ---- Plots ----

message("Generating QC plots...")
pdf(opt$out_plot, width = 12, height = 10)

tryCatch(
  print(CNAqc::plot_segments(x)),
  error = function(e) {
    message(paste("plot_segments failed:", conditionMessage(e)))
  }
)

tryCatch(
  print(CNAqc::plot_data_histogram(x)),
  error = function(e) {
    message(paste(
      "plot_data_histogram failed:", conditionMessage(e)
    ))
  }
)

if (opt$run_peak_analysis) {
  tryCatch(
    print(CNAqc::plot_peaks_analysis(x)),
    error = function(e) {
      message(paste(
        "plot_peaks_analysis failed:", conditionMessage(e)
      ))
    }
  )
}

dev.off()

# ---- Write metrics ----

message("Writing metrics TSV...")
write.table(
  score_df, opt$out_metrics,
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Done.")
