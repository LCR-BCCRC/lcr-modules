#!/usr/bin/env Rscript

# CNAqc analysis script
# Inputs: Battenberg subclones (filled), cellularity_ploidy, and MAF file
# Outputs: QC plots (PDF) and metrics/score TSV

suppressPackageStartupMessages({
    library(optparse)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--subclones"), type = "character", help = "Battenberg subclones.txt file (filled)"),
    make_option(c("--cellularity_ploidy"), type = "character", help = "Battenberg cellularity_ploidy.txt file"),
    make_option(c("--maf"), type = "character", help = "MAF file with somatic SNVs"),
    make_option(c("--tumour_id"), type = "character", help = "Tumour sample ID"),
    make_option(c("--out_plot"), type = "character", help = "Output PDF file for QC plots"),
    make_option(c("--out_metrics"), type = "character", help = "Output TSV file for QC metrics"),
    make_option(c("--min_depth"), type = "integer", default = 10, help = "Minimum read depth for variants [default: %default]"),
    make_option(c("--min_muts_per_segment"), type = "integer", default = 10, help = "Minimum mutations per segment for peak analysis [default: %default]"),
    make_option(c("--run_peak_analysis"), type = "logical", default = TRUE, help = "Run peak analysis [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
    description = "Run CNAqc QC analysis using Battenberg CNA and MAF SNVs")
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("subclones", "cellularity_ploidy", "maf", "tumour_id", "out_plot", "out_metrics")
for (arg in required_args) {
    if (is.null(opt[[arg]])) {
        stop(paste("Required argument missing:", arg))
    }
}

# Install CNAqc from GitHub if not already installed
if (!requireNamespace("CNAqc", quietly = TRUE)) {
    message("Installing CNAqc from GitHub (caravagnalab/CNAqc)...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github("caravagnalab/CNAqc", dependencies = TRUE, upgrade = "never")
}

suppressPackageStartupMessages({
    library(CNAqc)
    library(dplyr)
    library(ggplot2)
    library(data.table)
})

message("Loading input files...")

# ---- Load cellularity/ploidy ----
cp <- read.table(opt$cellularity_ploidy, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
purity <- as.numeric(cp$cellularity[1])
ploidy <- as.numeric(cp$ploidy[1])
message(sprintf("Purity: %.4f, Ploidy: %.4f", purity, ploidy))

# ---- Load Battenberg subclones ----
# Battenberg subclones format:
# chr  startpos  endpos  BAF  pval  logR  ntot  nMaj1_A  nMin1_A  frac1_A  [nMaj2_A  nMin2_A  frac2_A]  SDfrac_A  SDfrac_A_BS  SDlogR  [clonal fields]
sub <- fread(opt$subclones, sep = "\t", data.table = FALSE)

# Normalise column names (Battenberg uses varied capitalisations)
# Expected columns: chr, startpos, endpos, nMaj1_A, nMin1_A, frac1_A, [nMaj2_A, nMin2_A, frac2_A]
colnames(sub) <- tolower(colnames(sub))

# Check for required columns
required_sub_cols <- c("chr", "startpos", "endpos", "nmaj1_a", "nmin1_a", "frac1_a")
missing_cols <- required_sub_cols[!required_sub_cols %in% colnames(sub)]
if (length(missing_cols) > 0) {
    stop(paste("Battenberg subclones file missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Build CNA data.frame for CNAqc
# CNAqc expects columns: chr, from, to, Major, minor (integer copy numbers), CCF (optional)
# For subclonal segments, use the clone 1 (major clone) values
cna <- data.frame(
    chr = sub$chr,
    from = as.integer(sub$startpos),
    to = as.integer(sub$endpos),
    Major = as.integer(sub$nmaj1_a),
    minor = as.integer(sub$nmin1_a),
    CCF = as.numeric(sub$frac1_a),
    stringsAsFactors = FALSE
)

# Add chr prefix if not present (CNAqc expects chr1, chr2, ...)
if (!all(grepl("^chr", cna$chr))) {
    cna$chr <- paste0("chr", cna$chr)
}

# Drop any rows with NA copy numbers or zero-length segments
cna <- cna[!is.na(cna$Major) & !is.na(cna$minor) & (cna$to - cna$from) > 1, ]

message(sprintf("Loaded %d CNA segments from subclones file.", nrow(cna)))

# ---- Load MAF file ----
maf_raw <- fread(opt$maf, sep = "\t", skip = "#", data.table = FALSE)

# MAF standard column names (case-insensitive)
maf_cols <- tolower(colnames(maf_raw))
colnames(maf_raw) <- maf_cols

# Required columns for CNAqc SNVs: chr, from, to, ref, alt, DP, NV, VAF
# These map to: Chromosome, Start_Position, End_Position, Reference_Allele,
#               Tumor_Seq_Allele2, t_depth, t_alt_count
col_map <- list(
    chr = c("chromosome", "chrom"),
    from = c("start_position", "start"),
    ref = c("reference_allele"),
    alt = c("tumor_seq_allele2", "alt_allele"),
    dp = c("t_depth"),
    nv = c("t_alt_count")
)

find_col <- function(df_cols, candidates) {
    match <- candidates[candidates %in% df_cols]
    if (length(match) == 0) return(NULL)
    match[1]
}

chr_col <- find_col(maf_cols, col_map$chr)
from_col <- find_col(maf_cols, col_map$from)
ref_col <- find_col(maf_cols, col_map$ref)
alt_col <- find_col(maf_cols, col_map$alt)
dp_col <- find_col(maf_cols, col_map$dp)
nv_col <- find_col(maf_cols, col_map$nv)

for (col_name in c("chr_col", "from_col", "ref_col", "alt_col", "dp_col", "nv_col")) {
    if (is.null(get(col_name))) {
        stop(paste("Cannot find required MAF column for:", col_name))
    }
}

snvs <- data.frame(
    chr = maf_raw[[chr_col]],
    from = as.integer(maf_raw[[from_col]]),
    to = as.integer(maf_raw[[from_col]]) + 1L,
    ref = maf_raw[[ref_col]],
    alt = maf_raw[[alt_col]],
    DP = as.integer(maf_raw[[dp_col]]),
    NV = as.integer(maf_raw[[nv_col]]),
    stringsAsFactors = FALSE
)

# Compute VAF
snvs$VAF <- snvs$NV / snvs$DP

# Filter: SNVs only (single nucleotide substitutions), minimum depth
snvs <- snvs[
    !is.na(snvs$DP) & snvs$DP >= opt$min_depth &
    nchar(snvs$ref) == 1 & nchar(snvs$alt) == 1, ]

# Add chr prefix if not present
if (!all(grepl("^chr", snvs$chr))) {
    snvs$chr <- paste0("chr", snvs$chr)
}

message(sprintf("Loaded %d SNVs from MAF (after depth filter >= %d).", nrow(snvs), opt$min_depth))

# ---- Initialize CNAqc object ----
message("Initialising CNAqc object...")
x <- CNAqc::init(
    mutations = snvs,
    cna = cna,
    purity = purity,
    ref = "GRCh38"    # CNAqc handles grch37/hg38; adjust as needed
)

# ---- Peak analysis ----
if (opt$run_peak_analysis && nrow(snvs) >= opt$min_muts_per_segment) {
    message("Running peak analysis...")
    tryCatch({
        x <- CNAqc::analyze_peaks(x, matching_strategy = "closest")
    }, error = function(e) {
        message(sprintf("Peak analysis failed: %s", conditionMessage(e)))
    })
}

# ---- CNAqc score ----
message("Computing CNAqc score...")
score_df <- tryCatch({
    CNAqc::CNAqc_score(x)
}, error = function(e) {
    message(sprintf("CNAqc_score failed: %s", conditionMessage(e)))
    data.frame(score = NA_real_, QC = NA_character_, stringsAsFactors = FALSE)
})

if (!is.data.frame(score_df)) {
    score_df <- data.frame(score = as.numeric(score_df), QC = NA_character_, stringsAsFactors = FALSE)
}

score_df$tumour_id <- opt$tumour_id
score_df$purity    <- purity
score_df$ploidy    <- ploidy
score_df$n_snvs    <- nrow(snvs)
score_df$n_cna_segments <- nrow(cna)

message(sprintf("CNAqc score: %s", paste(capture.output(print(score_df)), collapse = " ")))

# ---- Plots ----
message("Generating QC plots...")
pdf(opt$out_plot, width = 12, height = 10)
tryCatch({
    p_cna <- CNAqc::plot_CNA(x)
    print(p_cna)
}, error = function(e) message(paste("plot_CNA failed:", conditionMessage(e))))

tryCatch({
    p_snvs <- CNAqc::plot_data_histogram(x)
    print(p_snvs)
}, error = function(e) message(paste("plot_data_histogram failed:", conditionMessage(e))))

if (opt$run_peak_analysis) {
    tryCatch({
        p_peaks <- CNAqc::plot_peaks_analysis(x)
        print(p_peaks)
    }, error = function(e) message(paste("plot_peaks_analysis failed:", conditionMessage(e))))
}
dev.off()

# ---- Write metrics ----
message("Writing metrics TSV...")
write.table(score_df, opt$out_metrics, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done.")
