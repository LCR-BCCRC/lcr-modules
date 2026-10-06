#!/usr/bin/env Rscript
# Builds a cohort-wide (not per-patient), per-locus lookup table of which known IPD-IMGT/HLA class
# I (A/B/C) alleles carry which amino acid at which position -- the reference data
# src/flag_hla_genotype_errors.R needs to check whether an apparent somatic substitution actually
# reconstructs a known germline allele state.
#
# Module-owned code. Not an upstream MHC Hammer script -- this module's own DNA-analysis arm never
# needed a multi-allele protein comparison; this exists purely for the opt-in genotype-error-
# detection feature (see options.hla_genotype_qc).
#
# Input is mhc_references/all_allele_info.csv -- confirmed (by directly downloading and inspecting
# the real Zenodo zip this module already fetches for typing/reference construction) to already
# contain, for every catalogued allele: allele_name (the ORIGINAL IMGT dotted form, e.g.
# "A*01:01:01:01" -- NOT this module's own lowercase/underscore contig-naming convention),
# gene (e.g. "HLA-A"), partial (boolean), and translation (the FULL PRECURSOR protein sequence,
# i.e. including the signal peptide, same numbering convention VEP itself uses since the per-
# allele GTF VEP annotates against has no separate signal-peptide feature).
#
# Real finding that shaped this script's design: a naive "exclude any allele whose translation
# length differs from the locus's modal length" approach (to justify plain same-length Hamming
# distance) was checked against the real data and found to exclude 4-7% of non-partial alleles per
# locus (hundreds of real, catalogued alleles) -- far too costly for a feature specifically meant
# to catch the rare correct-but-unassigned allele. Instead, EVERY non-partial allele is pairwise-
# aligned (Biostrings::pairwiseAlignment, global) against one fixed, data-driven per-locus
# reference sequence, so indel-bearing alleles are included too, just compared via the shared
# reference-anchored coordinate rather than raw sequence index.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument('--all_allele_info', nargs = 1, required = TRUE,
                    help = 'Path to mhc_references/all_allele_info.csv')
parser$add_argument('--translation_matrix_out', nargs = 1, required = TRUE,
                    help = 'Path to save the .rds lookup table (per-gene matrices + position maps)')
parser$add_argument('--position_states_out', nargs = 1, required = TRUE,
                    help = 'Path to save the long-format .csv.gz (gene, ref_position, amino_acid, group, allele_name)')
parser$add_argument('--excluded_alleles_out', nargs = 1, required = TRUE,
                    help = 'Path to save a CSV logging any allele excluded for a low alignment score')
parser$add_argument('--min_alignment_score_fraction', nargs = 1, type = 'double', default = 0.5,
                    help = 'An allele whose pairwiseAlignment score, divided by its own self-alignment score, falls below this fraction is excluded as likely-corrupt catalog data (default 0.5)')

args <- parser$parse_args()

GENES <- c("HLA-A", "HLA-B", "HLA-C")

cat("Reading", args$all_allele_info, "\n")
info <- fread(args$all_allele_info)
info <- info[gene %in% GENES & toupper(trimws(as.character(partial))) == "FALSE"]
cat("Rows after gene/partial filtering:", nrow(info), "\n")

data(BLOSUM62) # ships with Biostrings

# Allele-name group: everything before the first ":" (e.g. "A*03:01:01:01" -> "A*03"). Applied
# directly to the original IMGT dotted allele_name -- no lookup needed.
allele_group <- function(allele_name) sub(":.*$", "", allele_name)

translation_matrix <- list()
excluded_rows <- list()
position_states_rows <- list()

for (g in GENES) {

  sub_info <- info[gene == g]
  n_total <- nrow(sub_info)
  if (n_total == 0) {
    cat("WARNING: no non-partial alleles found for", g, "-- skipping\n")
    next
  }

  # Pick a single, data-driven reference-anchor sequence: among alleles at the modal translation
  # length, the most common EXACT sequence. Purely data-driven (no hardcoded allele names) so this
  # is robust across IMGT releases.
  lens <- nchar(sub_info$translation)
  modal_len <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
  modal_seqs <- sub_info$translation[lens == modal_len]
  ref_seq_str <- names(sort(table(modal_seqs), decreasing = TRUE))[1]
  ref_AA <- AAString(ref_seq_str)
  L_ref <- nchar(ref_seq_str)
  cat(sprintf("%s: %d non-partial alleles; modal length %d (n=%d); reference sequence length %d\n",
              g, n_total, modal_len, sum(lens == modal_len), L_ref))

  # Self-alignment score, used as the denominator for the exclusion sanity check below.
  self_aln <- pairwiseAlignment(pattern = ref_AA, subject = ref_AA, type = "global",
                                 substitutionMatrix = BLOSUM62, gapOpening = 10, gapExtension = 0.5)
  self_score <- score(self_aln)

  gene_letter <- sub("^HLA-", "", g)
  allele_mat <- matrix(NA_character_, nrow = nrow(sub_info), ncol = L_ref,
                        dimnames = list(sub_info$allele_name, NULL))
  own_to_ref_maps <- vector("list", nrow(sub_info))
  names(own_to_ref_maps) <- sub_info$allele_name

  n_excluded <- 0
  for (i in seq_len(nrow(sub_info))) {
    allele_name <- sub_info$allele_name[i]
    allele_AA <- AAString(sub_info$translation[i])

    aln <- pairwiseAlignment(pattern = allele_AA, subject = ref_AA, type = "global",
                              substitutionMatrix = BLOSUM62, gapOpening = 10, gapExtension = 0.5)

    if (score(aln) / self_score < args$min_alignment_score_fraction) {
      n_excluded <- n_excluded + 1
      excluded_rows[[length(excluded_rows) + 1]] <- data.table(
        gene = g, allele_name = allele_name, reason = "low_alignment_score",
        alignment_score = score(aln), self_score = self_score
      )
      next
    }

    aligned_pattern <- strsplit(as.character(alignedPattern(aln)), "")[[1]]
    aligned_subject <- strsplit(as.character(alignedSubject(aln)), "")[[1]]

    ref_pos <- 0L
    own_pos <- 0L
    own_to_ref <- integer(nchar(sub_info$translation[i]))
    for (col in seq_along(aligned_subject)) {
      subj_char <- aligned_subject[col]
      pat_char <- aligned_pattern[col]
      if (subj_char != "-") {
        ref_pos <- ref_pos + 1L
        allele_mat[i, ref_pos] <- pat_char # "-" recorded as a real gap state, not NA
      }
      if (pat_char != "-") {
        own_pos <- own_pos + 1L
        # If the reference itself has a gap here (this allele has an insertion relative to the
        # reference), there is no shared reference coordinate for this residue -- record NA.
        own_to_ref[own_pos] <- if (subj_char != "-") ref_pos else NA_integer_
      }
    }
    own_to_ref_maps[[allele_name]] <- own_to_ref
  }

  if (n_excluded > 0) {
    cat(sprintf("%s: excluded %d/%d alleles for a low alignment score (likely corrupt catalog data)\n",
                g, n_excluded, n_total))
  }

  translation_matrix[[gene_letter]] <- list(
    matrix = allele_mat,
    own_position_to_ref_position = own_to_ref_maps,
    reference_sequence = ref_seq_str,
    reference_length = L_ref
  )

  groups <- allele_group(rownames(allele_mat))
  for (ref_position in seq_len(L_ref)) {
    aa_col <- allele_mat[, ref_position]
    non_na <- !is.na(aa_col)
    if (!any(non_na)) next
    position_states_rows[[length(position_states_rows) + 1]] <- data.table(
      gene = gene_letter,
      ref_position = ref_position,
      amino_acid = aa_col[non_na],
      group = groups[non_na],
      allele_name = rownames(allele_mat)[non_na]
    )
  }
}

saveRDS(translation_matrix, args$translation_matrix_out)

position_states <- rbindlist(position_states_rows, use.names = TRUE)
fwrite(position_states, args$position_states_out)

if (length(excluded_rows) > 0) {
  fwrite(rbindlist(excluded_rows, use.names = TRUE), args$excluded_alleles_out)
} else {
  fwrite(data.table(gene = character(), allele_name = character(), reason = character(),
                     alignment_score = numeric(), self_score = numeric()),
         args$excluded_alleles_out)
}

cat("Wrote", args$translation_matrix_out, "\n")
cat("Wrote", args$position_states_out, "(", nrow(position_states), "rows)\n")
cat("Wrote", args$excluded_alleles_out, "\n")
