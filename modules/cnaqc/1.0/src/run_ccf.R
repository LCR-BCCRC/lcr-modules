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

CCF_COLS <- c("karyotype", "segment_id", "QC_PASS",
              "mutation_multiplicity", "CCF_estimate",
              "CCF_CI_low", "CCF_CI_high")

BT_COLS <- c("bt_major_1", "bt_minor_1", "bt_frac_1", "bt_karyotype_1",
             "bt_major_2", "bt_minor_2", "bt_frac_2", "bt_karyotype_2")

write_all_na <- function(path, tumour_id, muts) {
  for (col in c(CCF_COLS, BT_COLS)) {
    if (!col %in% colnames(muts)) muts[[col]] <- NA
  }
  out <- cbind(tumour_id = tumour_id, muts)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
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

snvs$VAF  <- snvs$NV / snvs$DP
snvs$type <- ifelse(nchar(snvs$ref) == 1L & nchar(snvs$alt) == 1L, "SNV", "indel")
snvs <- snvs[!is.na(snvs$DP) & snvs$DP >= opt$min_depth, ]

if (!all(grepl("^chr", snvs$chr))) snvs$chr <- paste0("chr", snvs$chr)
message(sprintf("Loaded %d mutations (depth >= %d): %d SNVs, %d indels.",
  nrow(snvs), opt$min_depth,
  sum(snvs$type == "SNV"), sum(snvs$type == "indel")))

# ---- CNAqc init ----

x <- CNAqc::init(mutations = snvs, cna = cna, purity = purity, ref = opt$ref)

# ---- Peak analysis (required before compute_CCF) ----

if (nrow(snvs) < opt$min_muts_per_segment) {
  message("Too few mutations for peak analysis — writing all mutations with NA CCF.")
  write_all_na(opt$out_ccf, opt$tumour_id, snvs)
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
  message("Writing all mutations with NA CCF.")
  write_all_na(opt$out_ccf, opt$tumour_id, snvs)
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
  message("Writing all mutations with NA CCF.")
  write_all_na(opt$out_ccf, opt$tumour_id, snvs)
  quit(save = "no", status = 0)
}

# ---- Export (all input mutations, with NA for those CNAqc could not process) ----

join_key <- c("chr", "from", "ref", "alt")

# All clonal mutations (simple and complex karyotypes)
clonal_all <- tryCatch(CNAqc::Mutations(x, cna = "clonal"), error = function(e) NULL)

# CCF estimates — only available for karyotypes where CNAqc fitted peaks
ccf_estimates <- tryCatch(CNAqc::CCF(x), error = function(e) {
  message(paste("CCF() getter failed:", conditionMessage(e)))
  NULL
})

# Left-join CCF estimates onto all clonal mutations so high-CN segments
# (e.g. 3:1, 5:0, 11:0) get karyotype/segment_id but NA for CCF/multiplicity.
# NULL (CCF() errored) and empty (no estimates produced) are handled separately:
# when empty we still know the column schema and can add NA columns to keep
# the output schema consistent with the success path.
if (!is.null(clonal_all) && nrow(clonal_all) > 0 && !is.null(ccf_estimates)) {
  ccf_cols <- setdiff(colnames(ccf_estimates), colnames(clonal_all))
  if (nrow(ccf_estimates) > 0) {
    clonal_df <- dplyr::left_join(
      clonal_all,
      ccf_estimates[, c(join_key, ccf_cols), drop = FALSE],
      by = join_key
    )
  } else {
    clonal_df <- clonal_all
    for (col in ccf_cols) clonal_df[[col]] <- NA
  }
} else {
  clonal_df <- clonal_all
}

# Group 2: subclonal mutations — mapped to a segment but no CCF estimate
subclonal_df <- tryCatch(CNAqc::Mutations(x, cna = "subclonal"), error = function(e) NULL)

# Combine all mutations CNAqc mapped to a segment
mapped_df <- dplyr::bind_rows(clonal_df, subclonal_df)

# Group 3: unmapped mutations — not overlapping any CNA segment
if (!is.null(mapped_df) && nrow(mapped_df) > 0) {
  unmapped_df <- dplyr::anti_join(snvs, mapped_df, by = join_key)
} else {
  unmapped_df <- snvs
}

# Fill NA for all CNAqc-derived columns in unmapped mutations
if (nrow(unmapped_df) > 0 && !is.null(mapped_df) && nrow(mapped_df) > 0) {
  for (col in setdiff(colnames(mapped_df), colnames(unmapped_df))) {
    unmapped_df[[col]] <- NA
  }
}

out_df <- dplyr::bind_rows(mapped_df, unmapped_df)

# Guarantee CCF-derived columns are always present, even when CCF() errored
# and we couldn't read the schema from ccf_estimates.
for (col in c(CCF_COLS, BT_COLS)) {
  if (!col %in% colnames(out_df)) out_df[[col]] <- NA
}

# ---- Re-join Battenberg segment CN data ----
# CNAqc leaves karyotype=NA for subclonal mutations (the segment is a CN mixture
# so no single karyotype can be assigned). Re-joining to the original Battenberg
# data recovers both clones' CN values for all segment-mapped mutations.

# Ensure clone 2 columns exist (absent when Battenberg found no subclonality)
for (col in c("nmaj2_a", "nmin2_a", "frac2_a")) {
  if (!col %in% colnames(sub)) sub[[col]] <- NA
}

sub_cn_annot <- data.frame(
  .chr           = if_else(grepl("^chr", sub$chr), sub$chr, paste0("chr", sub$chr)),
  .startpos      = as.integer(sub$startpos),
  .endpos        = as.integer(sub$endpos),
  bt_major_1     = as.integer(sub$nmaj1_a),
  bt_minor_1     = as.integer(sub$nmin1_a),
  bt_frac_1      = as.numeric(sub$frac1_a),
  bt_karyotype_1 = paste0(sub$nmaj1_a, ":", sub$nmin1_a),
  bt_major_2     = as.integer(sub$nmaj2_a),
  bt_minor_2     = as.integer(sub$nmin2_a),
  bt_frac_2      = as.numeric(sub$frac2_a),
  bt_karyotype_2 = ifelse(
    !is.na(sub$nmaj2_a) & !is.na(sub$nmin2_a),
    paste0(sub$nmaj2_a, ":", sub$nmin2_a),
    NA_character_
  ),
  stringsAsFactors = FALSE
)

# Parse chr:start:end from segment_id (format: chr:start:end:Major:minor:CCF)
out_df <- out_df %>%
  mutate(
    .seg_chr   = if_else(!is.na(segment_id), sub(":.*", "", segment_id), NA_character_),
    .seg_start = if_else(!is.na(segment_id),
                         suppressWarnings(as.integer(gsub("^[^:]+:([^:]+):.*", "\\1", segment_id))),
                         NA_integer_),
    .seg_end   = if_else(!is.na(segment_id),
                         suppressWarnings(as.integer(gsub("^[^:]+:[^:]+:([^:]+):.*", "\\1", segment_id))),
                         NA_integer_)
  ) %>%
  dplyr::left_join(sub_cn_annot,
                   by = c(".seg_chr"   = ".chr",
                          ".seg_start" = ".startpos",
                          ".seg_end"   = ".endpos")) %>%
  select(-.seg_chr, -.seg_start, -.seg_end)

out_df <- cbind(tumour_id = opt$tumour_id, out_df)

message(sprintf(
  "Writing %d mutations: %d clonal with CCF, %d subclonal, %d unmapped.",
  nrow(out_df),
  if (!is.null(clonal_df)) nrow(clonal_df) else 0L,
  if (!is.null(subclonal_df)) nrow(subclonal_df) else 0L,
  nrow(unmapped_df)
))
write.table(out_df, opt$out_ccf, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done.")
