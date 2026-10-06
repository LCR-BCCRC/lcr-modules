#!/usr/bin/env Rscript
# Flags apparent somatic missense mutations in mhc_hammer's own {patient_id}_mutations.csv that
# look like HLA-HD typing-assignment errors rather than real somatic mutations: if HLA-HD assigns
# the wrong (but closely related) germline allele, every position where that wrong allele differs
# from the patient's TRUE allele shows up as a spurious "somatic" substitution. This script never
# changes or drops a mutation call -- it only adds a parallel, purely informational annotation file
# for a human to review. See options.hla_genotype_qc in default.yaml for the full rationale and the
# real HLA-A p.K292E / mature K268E prototype case this was built to catch.
#
# Module-owned code. Works identically for either of this module's two mutation-calling arms
# (paired DNA, or tumour-only) since their own {patient_id}_mutations.csv schemas are identical by
# design -- the caller just passes --pathway for provenance.
#
# Core principle (the user's own framing): a putative HLA somatic mutation is suspicious whenever
# the observed "mutant" state reconstructs a KNOWN GERMLINE HLA allele, especially when several
# apparent mutations on the SAME assigned allele collectively reconstruct the same alternate
# allele. Every apparent mutation on a given (tumour_sample_name, assigned allele) is evaluated
# jointly against the full candidate pool, not independently.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument('--mutations_csv', nargs = 1, required = TRUE,
                    help = 'Either arm\'s {patient_id}_mutations.csv -- identical schema')
parser$add_argument('--translation_matrix', nargs = 1, required = TRUE,
                    help = 'Path to the .rds built by build_allele_lookup.R')
parser$add_argument('--sample_genotype', nargs = '+', required = TRUE,
                    help = 'One or more tumour_sample_name=path_to_hla_alleles.csv tokens')
parser$add_argument('--pathway', nargs = 1, required = TRUE, choices = c('paired', 'tumour_only'))
parser$add_argument('--imgt_release', nargs = 1, required = TRUE)
parser$add_argument('--signal_peptide_length_a', nargs = 1, type = 'integer', default = 24)
parser$add_argument('--signal_peptide_length_b', nargs = 1, type = 'integer', default = 24)
parser$add_argument('--signal_peptide_length_c', nargs = 1, type = 'integer', default = 24)
parser$add_argument('--max_ref_mismatch_fraction', nargs = 1, type = 'double', default = 0.5,
                    help = 'Abort if the fraction of rows where VEP\'s own REF disagrees with the catalog exceeds this')
parser$add_argument('--mutation_save_path', nargs = 1, required = TRUE)

args <- parser$parse_args()

SIGNAL_PEPTIDE_LENGTH <- c(A = args$signal_peptide_length_a, B = args$signal_peptide_length_b,
                           C = args$signal_peptide_length_c)

# --- parse --sample_genotype tokens into sample_id -> hla_alleles.csv path, then load each
# distinct genotype file once (the paired arm passes the same patient-level file repeated per
# tumour sample; the tumour-only arm passes one independently-typed file per sample -- confirmed
# in this module's own README, "Per-tumour typing, not shared per-patient") -----------------------
genotype_tokens <- strsplit(args$sample_genotype, "=", fixed = TRUE)
sample_to_genotype_path <- setNames(
  vapply(genotype_tokens, `[[`, character(1), 2),
  vapply(genotype_tokens, `[[`, character(1), 1)
)
genotype_cache <- new.env()
load_genotype <- function(path) {
  if (is.null(genotype_cache[[path]])) {
    dt <- fread(path)
    if (!all(c("gene", "allele1", "allele2") %in% names(dt))) {
      setnames(dt, c("gene", "allele1", "allele2")[seq_len(ncol(dt))])
    }
    genotype_cache[[path]] <- dt
  }
  genotype_cache[[path]]
}

# hla_a_03_01_01_01 -> list(gene_letter = "A", dotted = "A*03:01:01:01") -- exact reverse of
# upstream's own contig-naming transform (create_mhc_hammer_references.R:125:
# paste0("hla_", tolower(gsub('\\*|:', '_', allele)))).
parse_contig_name <- function(chrom) {
  stripped <- sub("^hla_", "", chrom, ignore.case = TRUE)
  parts <- strsplit(stripped, "_", fixed = TRUE)[[1]]
  gene_letter <- toupper(parts[1])
  dotted <- paste0(gene_letter, "*", paste(parts[-1], collapse = ":"))
  list(gene_letter = gene_letter, dotted = dotted)
}

allele_group <- function(dotted_allele_name) sub(":.*$", "", dotted_allele_name)

# Gap-aware distance between two rows of the same (allele x ref_position) matrix: both gaps at a
# position is NOT a mismatch (both alleles simply lack that residue); one gap and one real residue
# IS a mismatch; two different real residues IS a mismatch. Written to avoid any NA propagation
# through `sum()`.
gapaware_distance <- function(vec_a, vec_b) {
  both_gap <- is.na(vec_a) & is.na(vec_b)
  one_gap <- xor(is.na(vec_a), is.na(vec_b))
  differ <- !is.na(vec_a) & !is.na(vec_b) & vec_a != vec_b
  sum(one_gap | differ)
}

cat("Loading translation matrix from", args$translation_matrix, "\n")
translation_matrix <- readRDS(args$translation_matrix)

cat("Loading", args$mutations_csv, "\n")
muts <- fread(args$mutations_csv)

qc_cols <- c("tumour_sample_name", "pathway", "locus", "assigned_allele", "vep_protein_position",
             "ipd_mature_position", "ref_matches_catalog", "vep_amino_acids",
             "other_patient_allele_matches_alt", "alt_known_germline_same_group",
             "alt_known_germline_alleles_same_group", "alt_known_germline_other_groups",
             "alt_known_germline_alleles_other_groups", "best_alternate_allele",
             "cross_group_candidate", "n_calls_explained_by_alternate", "n_calls_conflicting",
             "candidate_region_depth_adequate", "n_conflicting_observed_sites",
             "sequence_distance_to_assigned", "alt_fraction", "suspicion_level", "reason",
             "imgt_release")

if (nrow(muts) == 0) {
  fwrite(data.table(matrix(ncol = length(qc_cols), nrow = 0, dimnames = list(NULL, qc_cols))),
         args$mutation_save_path)
  cat("No mutation rows -- wrote an empty (header-only) QC file.\n")
  quit(save = "no", status = 0)
}

blank_row <- function(sample_name, reason) {
  row <- as.list(rep(NA, length(qc_cols)))
  names(row) <- qc_cols
  row$tumour_sample_name <- sample_name
  row$suspicion_level <- "none"
  row$reason <- reason
  row$pathway <- args$pathway
  row$imgt_release <- args$imgt_release
  row
}

is_simple_missense <- grepl("^[A-Z]/[A-Z]$", muts$vep_amino_acids) &
  grepl("^[0-9]+$", trimws(as.character(muts$vep_protein_position)))

out_rows <- vector("list", nrow(muts))
ref_positions <- rep(NA_integer_, nrow(muts)) # kept alongside qc for the joint-ranking pass below

for (i in seq_len(nrow(muts))) {
  if (!is_simple_missense[i]) {
    out_rows[[i]] <- blank_row(muts$tumour_sample_name[i],
                               "not a simple missense substitution -- out of scope for this check")
    next
  }

  contig <- parse_contig_name(muts$chrom[i])
  gene_letter <- contig$gene_letter
  assigned_allele <- contig$dotted
  gene_lookup <- translation_matrix[[gene_letter]]

  if (is.null(gene_lookup) || !(assigned_allele %in% rownames(gene_lookup$matrix))) {
    r <- blank_row(muts$tumour_sample_name[i], sprintf(
      "assigned allele %s not found in the reference lookup table (excluded for a low alignment score, or an IMGT release mismatch -- see excluded_alleles.csv)",
      assigned_allele
    ))
    r$locus <- gene_letter
    r$assigned_allele <- assigned_allele
    out_rows[[i]] <- r
    next
  }

  own_to_ref <- gene_lookup$own_position_to_ref_position[[assigned_allele]]
  vep_pos <- as.integer(muts$vep_protein_position[i])
  ref_position <- if (!is.na(vep_pos) && vep_pos >= 1 && vep_pos <= length(own_to_ref)) {
    own_to_ref[vep_pos]
  } else {
    NA_integer_
  }

  r <- blank_row(muts$tumour_sample_name[i], NA_character_)
  r$locus <- gene_letter
  r$assigned_allele <- assigned_allele
  r$vep_protein_position <- vep_pos

  if (is.na(ref_position)) {
    r$reason <- "assigned allele's own position maps to an insertion with no shared reference coordinate -- cannot be compared"
    out_rows[[i]] <- r
    next
  }
  ref_positions[i] <- ref_position
  r$ipd_mature_position <- ref_position - SIGNAL_PEPTIDE_LENGTH[[gene_letter]]

  aa_parts <- strsplit(muts$vep_amino_acids[i], "/", fixed = TRUE)[[1]]
  ref_aa <- aa_parts[1]
  alt_aa <- aa_parts[2]
  r$vep_amino_acids <- muts$vep_amino_acids[i]

  # unname() matters here: matrix[rowname, col_index] carries the rowname over as a `names`
  # attribute on the scalar it returns (R keeps whichever dimnames survive the drop), which makes
  # identical()/== against a plain string look like a mismatch unless stripped first.
  catalog_ref_aa <- unname(gene_lookup$matrix[assigned_allele, ref_position])
  r$ref_matches_catalog <- !is.na(catalog_ref_aa) && catalog_ref_aa == ref_aa

  # --- germline-state check: every allele carrying ALT at this ref_position ----------------------
  col <- gene_lookup$matrix[, ref_position]
  carriers <- names(col)[!is.na(col) & col == alt_aa]
  assigned_group <- allele_group(assigned_allele)
  carrier_groups <- allele_group(carriers)
  same_group_carriers <- carriers[carrier_groups == assigned_group]
  other_group_carriers <- carriers[carrier_groups != assigned_group]

  r$alt_known_germline_same_group <- length(same_group_carriers) > 0
  r$alt_known_germline_alleles_same_group <- paste(same_group_carriers, collapse = ";")
  r$alt_known_germline_other_groups <- length(other_group_carriers) > 0
  r$alt_known_germline_alleles_other_groups <- paste(other_group_carriers, collapse = ";")

  # --- patient's other allele at this locus, cross-check ------------------------------------------
  genotype_path <- sample_to_genotype_path[[muts$tumour_sample_name[i]]]
  other_match <- NA
  if (!is.null(genotype_path)) {
    genotype <- load_genotype(genotype_path)
    g_row <- genotype[toupper(gsub("^HLA-", "", gene)) == gene_letter]
    if (nrow(g_row) >= 1) {
      typed_alleles <- c(g_row$allele1[1], g_row$allele2[1])
      typed_alleles <- typed_alleles[!is.na(typed_alleles) & nzchar(typed_alleles)]
      other_allele <- typed_alleles[typed_alleles != assigned_allele]
      if (length(other_allele) >= 1) {
        other_allele <- other_allele[1]
        if (other_allele %in% rownames(gene_lookup$matrix)) {
          other_aa <- gene_lookup$matrix[other_allele, ref_position]
          other_match <- !is.na(other_aa) && other_aa == alt_aa
        }
      }
    }
  }
  r$other_patient_allele_matches_alt <- other_match

  tumour_alt_dp <- suppressWarnings(as.numeric(muts$tumour_alt_dp[i]))
  tumour_ref_dp <- suppressWarnings(as.numeric(muts$tumour_ref_dp[i]))
  r$alt_fraction <- if (!is.na(tumour_alt_dp) && !is.na(tumour_ref_dp) && (tumour_alt_dp + tumour_ref_dp) > 0) {
    tumour_alt_dp / (tumour_alt_dp + tumour_ref_dp)
  } else {
    NA_real_
  }

  # --- coarse depth-adequacy proxy (per the resolved v1 decision -- precise per-site pileup
  # conflict counting is deferred) ------------------------------------------------------------------
  tumour_dp <- suppressWarnings(as.numeric(muts$tumour_dp[i]))
  r$candidate_region_depth_adequate <- !is.na(tumour_dp) && tumour_dp > 0
  r$n_conflicting_observed_sites <- NA_integer_

  out_rows[[i]] <- r
}

qc <- rbindlist(out_rows, use.names = TRUE, fill = TRUE)
qc[, ref_position := ref_positions]

# These four columns are never set in the per-row loop above -- only later, in the joint-ranking
# `:=` pass -- so every row's value at this point is still the generic logical NA blank_row() gave
# it, and rbindlist() infers a plain `logical` column type from that. Pin down the real types here
# so the later assignment doesn't force data.table to silently coerce (and corrupt) real
# character/integer/numeric values into a logical column. `:=` updates an EXISTING column in
# place and preserves its current type (coercing the RHS to match, not the other way round) --
# wholesale `$<-` replacement is what actually changes a column's type.
qc$best_alternate_allele <- rep(NA_character_, nrow(qc))
qc$n_calls_explained_by_alternate <- rep(NA_integer_, nrow(qc))
qc$n_calls_conflicting <- rep(NA_integer_, nrow(qc))
qc$sequence_distance_to_assigned <- rep(NA_real_, nrow(qc))

# --- self-consistency abort check: a few isolated ref_matches_catalog == FALSE rows are plausibly
# another instance of the exact phenomenon this script exists to catch; a majority-mismatching
# patient indicates a real numbering/parsing bug that needs fixing first, not data to report ------
checkable <- qc[!is.na(ref_matches_catalog)]
if (nrow(checkable) > 0) {
  mismatch_fraction <- sum(!checkable$ref_matches_catalog) / nrow(checkable)
  cat(sprintf("ref_matches_catalog: %d/%d rows mismatch (%.1f%%)\n",
              sum(!checkable$ref_matches_catalog), nrow(checkable), 100 * mismatch_fraction))
  if (mismatch_fraction > args$max_ref_mismatch_fraction) {
    stop(sprintf(
      "%.1f%% of checkable rows have a REF amino acid that doesn't match the assigned allele's own catalog translation (threshold: %.1f%%) -- this is far more consistent with a numbering/parsing bug than with real biology. Aborting rather than reporting unreliable QC output. See --max_ref_mismatch_fraction to adjust.",
      100 * mismatch_fraction, 100 * args$max_ref_mismatch_fraction
    ))
  }
}

# --- joint candidate ranking per (tumour_sample_name, locus, assigned_allele) group, using ALL of
# that group's apparent mutations together -- the core of the user's own design -------------------
qc[, row_idx := .I]
qc[, group_key := paste(tumour_sample_name, locus, assigned_allele, sep = "||")]
# Every row with a resolved ref_position enters the loop below, even ones with no germline carrier
# anywhere at this position -- the in-loop `candidate_pool` check (not this filter) is what decides
# and records that outcome. An earlier version filtered rows out here when both carrier columns
# were empty, which meant a solo mutation with truly no known germline match at this locus never
# ran the loop body at all and was left with an unset (NA) `reason`.
comparable <- qc[!is.na(ref_position)]

for (key in unique(comparable$group_key)) {
  rows_idx <- qc[group_key == key & !is.na(ref_position), row_idx]
  if (length(rows_idx) == 0) next
  group_rows <- qc[rows_idx]

  candidate_pool <- unique(unlist(strsplit(
    group_rows$alt_known_germline_alleles_same_group[nzchar(group_rows$alt_known_germline_alleles_same_group)],
    ";", fixed = TRUE
  )))
  cross_group <- FALSE
  if (length(candidate_pool) == 0) {
    candidate_pool <- unique(unlist(strsplit(
      group_rows$alt_known_germline_alleles_other_groups[nzchar(group_rows$alt_known_germline_alleles_other_groups)],
      ";", fixed = TRUE
    )))
    cross_group <- length(candidate_pool) > 0
  }
  if (length(candidate_pool) == 0) {
    qc[rows_idx, `:=`(suspicion_level = "none",
                        reason = "ALT is not a known germline state anywhere at this locus")]
    next
  }

  gene_letter <- group_rows$locus[1]
  gene_lookup <- translation_matrix[[gene_letter]]
  assigned_allele <- group_rows$assigned_allele[1]
  assigned_vec <- gene_lookup$matrix[assigned_allele, ]

  row_rps <- qc$ref_position[rows_idx]
  row_alts <- vapply(qc$vep_amino_acids[rows_idx], function(x) strsplit(x, "/", fixed = TRUE)[[1]][2], character(1))
  row_refs <- vapply(qc$vep_amino_acids[rows_idx], function(x) strsplit(x, "/", fixed = TRUE)[[1]][1], character(1))

  scores <- data.table(candidate = candidate_pool)
  scores[, n_explained := vapply(candidate, function(cand) {
    cand_vec <- gene_lookup$matrix[cand, ]
    cand_aa <- cand_vec[row_rps]
    sum(!is.na(cand_aa) & cand_aa == row_alts)
  }, integer(1))]
  scores[, n_conflicting := vapply(candidate, function(cand) {
    cand_vec <- gene_lookup$matrix[cand, ]
    cand_aa <- cand_vec[row_rps]
    sum(!is.na(cand_aa) & cand_aa != row_alts & cand_aa != row_refs)
  }, integer(1))]
  scores[, distance := vapply(candidate, function(cand) {
    gapaware_distance(assigned_vec, gene_lookup$matrix[cand, ])
  }, numeric(1))]

  setorder(scores, -n_explained, n_conflicting, distance, candidate)
  best <- scores[1]

  qc[rows_idx, `:=`(
    best_alternate_allele = best$candidate,
    cross_group_candidate = cross_group,
    n_calls_explained_by_alternate = best$n_explained,
    n_calls_conflicting = best$n_conflicting,
    sequence_distance_to_assigned = best$distance,
    suspicion_level = fcase(
      best$n_explained >= 2 & best$n_conflicting == 0, "high",
      best$n_explained >= 1, "medium",
      default = "none"
    ),
    reason = sprintf(
      "ALT reconstructs %s (explains %d/%d apparent mutation(s) on this allele, %d conflicting site(s), sequence distance %d)%s",
      best$candidate, best$n_explained, length(rows_idx), best$n_conflicting, best$distance,
      if (cross_group) " [candidate is NOT in the same first-field group as the assigned allele]" else ""
    )
  )]
}

qc[, c("row_idx", "group_key", "ref_position") := NULL]
fwrite(qc[, ..qc_cols], args$mutation_save_path)
cat("Wrote", nrow(qc), "QC rows to", args$mutation_save_path, "\n")
