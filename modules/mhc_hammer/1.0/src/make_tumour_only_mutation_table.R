#!/usr/bin/env Rscript
# Builds a per-patient mutation table for mhc_hammer's tumour-only (no matched normal DNA)
# pathway, in a column-compatible shape with upstream MHC Hammer's own bin/make_mutation_table.R
# output -- so that src/mutations_to_maf.R (this module's own downstream script) can consume either
# one unmodified.
#
# This is module-owned code, NOT a copy or adaptation of upstream's bin/make_mutation_table.R
# (off-limits to redistribute/modify under its Cancer Research Horizons academic-use licence -- see
# the licensing note near the top of mhc_hammer.smk). It exists because that script cannot run in
# tumour-only mode at all, confirmed by reading it directly: `--wxs_gl_bam_files` is
# `required=TRUE`; it `stop("There should be a tumour and germline column")`s if a VEP table is
# missing either sample's AD/DP genotype column (which a tumour-only Mutect2 VCF never has, since
# Mutect2 only emits a genotype column for samples passed via `-I`); and it does a hard
# `merge()` + `stop("Why are some tumours missing matched germlines?")` against a
# `normal_sample_name` column from an `inventory.csv` this pathway has no equivalent of. This
# script reimplements only the non-germline-dependent parts of that logic (VEP consequence
# splitting, mutation-type derivation, per-mutation BAM read-count backfill via upstream's own
# mutation_table_function.R, called by -- not copied -- the same way every other rule in this file
# calls upstream's scripts), and fills every germline_* output column with NA so the output schema
# matches upstream's own column set exactly.
#
# Each tumour sample passed in was independently typed and referenced from its OWN reads (this
# module's tumour-only pathway does not collapse multiple tumours from one patient to a single
# consensus type) -- tumour_sample_name distinguishes rows so mutation calls can be compared across
# tumours/biopsies from the same patient, each evaluated against its own personalised reference.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(deepSNV))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument('--vep_tables',
                    help = 'Path to VEP txt files.',
                    nargs = "+",
                    required = TRUE)

parser$add_argument('--tumour_bam_files',
                    help = 'Path to tumour-only allele-filtered BAMs (across every tumour sample for this patient).',
                    nargs = "+",
                    required = TRUE)

parser$add_argument('--mutation_save_path',
                    help = 'Path to save mutation file as',
                    required = TRUE)

parser$add_argument('--scripts_dir', nargs = 1,
                    help = 'Path to mhc_hammer_scripts_dir/bin -- used only to source() upstream\'s own mutation_table_function.R for get_bam_read_count().',
                    required = TRUE)

args <- parser$parse_args()
vep_tables <- args$vep_tables
tumour_bam_files <- args$tumour_bam_files
mutation_save_path <- args$mutation_save_path
scripts_dir <- args$scripts_dir

source(paste0(scripts_dir, "/mutation_table_function.R"))

cat("vep_tables=", vep_tables, "\n")
cat("tumour_bam_files=", tumour_bam_files, "\n")
cat("mutation_save_path=", mutation_save_path, "\n")
cat("scripts_dir=", scripts_dir, "\n")

tumour_samples <- unique(gsub("_wxs_novoalign.*", "", tumour_bam_files))

# first get a unique set of all mutations (identical logic to upstream's own tumour-side handling,
# just without requiring a germline AD/DP column to also be present in each VEP table)
tumour_mutations <- data.table()
for (vep_tab in vep_tables) {

  region_mutations <- fread(vep_tab)

  tumour_dp_col <- names(region_mutations)[names(region_mutations) %in% paste0(tumour_samples, ".DP")]
  tumour_ad_col <- names(region_mutations)[names(region_mutations) %in% paste0(tumour_samples, ".AD")]

  if (length(tumour_dp_col) != 1 | length(tumour_ad_col) != 1) {
    stop("There should be exactly one tumour DP/AD column pair in a tumour-only VEP table")
  }

  setnames(region_mutations,
           old = c(tumour_dp_col, tumour_ad_col),
           new = c("tumour_dp", "tumour_ad"))

  region_mutations[, tumour_sample_name := gsub(".DP$", "", tumour_dp_col)]

  tumour_mutations <- rbindlist(list(tumour_mutations, region_mutations), use.names = TRUE)

}

tumour_mutations[, tumour_ref_dp := gsub(",.*", "", tumour_ad)]
tumour_mutations[, tumour_alt_dp := gsub(".*,", "", tumour_ad)]

# No germline sample exists for this pathway -- every germline_* field is NA by construction, not
# derived from any column in the VEP table. germline_sample_name is set later, only on
# tumour_bam_read_counts (mirrors upstream's own structure, where that column only ever exists on
# the bam-read-count table, never on tumour_mutations -- setting it here too would collide at the
# merge below, since both tables would then carry the same column name).
tumour_mutations[, germline_dp := NA_character_]
tumour_mutations[, germline_ad := NA_character_]
tumour_mutations[, germline_ref_dp := NA_character_]
tumour_mutations[, germline_alt_dp := NA_character_]

# split the VEP consequence line into separate columns (identical logic to upstream)
for (line_idx in 1:nrow(tumour_mutations)) {
  vep_line <- tumour_mutations[line_idx]$CSQ
  vep_line_split <- strsplit(x = vep_line, split = "|", fixed = TRUE)[[1]]
  tumour_mutations[line_idx, vep_impact := vep_line_split[3]]
  tumour_mutations[line_idx, vep_feature_type := vep_line_split[4]]
  tumour_mutations[line_idx, vep_feature := vep_line_split[5]]
  tumour_mutations[line_idx, vep_exon := vep_line_split[6]]
  tumour_mutations[line_idx, vep_intron := vep_line_split[7]]
  tumour_mutations[line_idx, vep_cdna_position := vep_line_split[8]]
  tumour_mutations[line_idx, vep_cds_position := vep_line_split[9]]
  tumour_mutations[line_idx, vep_protein_position := vep_line_split[10]]
  tumour_mutations[line_idx, vep_amino_acids := vep_line_split[11]]
  tumour_mutations[line_idx, vep_codons := vep_line_split[12]]
  tumour_mutations[line_idx, vep_existing_variation := vep_line_split[13]]
  tumour_mutations[line_idx, vep_existing_distance := vep_line_split[14]]
  tumour_mutations[line_idx, vep_existing_strand := vep_line_split[15]]
  tumour_mutations[line_idx, vep_existing_flag := vep_line_split[16]]

  vep_consequence_line <- vep_line_split[2]
  vep_consequence_split <- strsplit(vep_consequence_line, "&")[[1]]
  vep_consequence_split <- sort(vep_consequence_split)
  new_vep_consequence <- paste0(vep_consequence_split, collapse = "&")
  tumour_mutations[line_idx, vep_consequence := new_vep_consequence]
}

# mutation type + synthetic start/end (identical logic to upstream)
tumour_mutations[, nchar_ref := nchar(REF)]
tumour_mutations[, nchar_alt := nchar(ALT)]
tumour_mutations[nchar_ref == 1 & nchar_alt == 1, mut_type := "SNV"]
tumour_mutations[nchar_ref > nchar_alt & nchar_alt == 1, mut_type := "DEL"]
tumour_mutations[nchar_ref < nchar_alt & nchar_ref == 1, mut_type := "INS"]
tumour_mutations[nchar_ref == nchar_alt & nchar_alt > 1, mut_type := "MNV"]
tumour_mutations[is.na(mut_type), mut_type := "COMPLEX"]

tumour_mutations[, start_dna := POS]
tumour_mutations[mut_type == "MNV", end_dna := POS + nchar_alt - 1]
tumour_mutations[mut_type == "DEL", end_dna := POS + nchar_ref - 1]

tumour_mutations[, c("nchar_ref", "nchar_alt") := NULL]

# loop over the tumour samples and get the ref/alt/N counts on the unique set of mutations
# (identical logic to upstream's own tumour-side loop -- no germline-side loop exists here)
unique_mutations <- unique(tumour_mutations[, c("CHROM", "POS", "REF", "ALT",
                                                "mut_type", "start_dna", "end_dna")])
tumour_bam_read_counts <- data.table()
for (tumour_sample in tumour_samples) {

  region_bam_read_counts <- copy(unique_mutations)
  region_bam_read_counts[, tumour_sample_name := tumour_sample]
  region_bam_read_counts[, bam_path := paste0(tumour_sample, "_wxs_novoalign.",
                                              CHROM, ".sorted.filtered.bam")]

  for (line_idx in 1:nrow(region_bam_read_counts)) {

    mutation_type <- region_bam_read_counts[line_idx]$mut_type
    chr <- region_bam_read_counts[line_idx]$CHROM
    start <- region_bam_read_counts[line_idx]$start_dna
    stop <- region_bam_read_counts[line_idx]$end_dna
    ref <- region_bam_read_counts[line_idx]$REF
    alt <- region_bam_read_counts[line_idx]$ALT

    if (mutation_type != "COMPLEX") {
      bam_read_count_output <- get_bam_read_count(bam_path = region_bam_read_counts[line_idx]$bam_path,
                                                  chr, start, stop, ref, alt, mutation_type)

      region_bam_read_counts[line_idx, tumour_ref_count := bam_read_count_output$bam_ref_count]
      region_bam_read_counts[line_idx, tumour_alt_count := bam_read_count_output$bam_alt_count]
      region_bam_read_counts[line_idx, tumour_N_count := bam_read_count_output$bam_N_count]
    } else {
      region_bam_read_counts[line_idx, tumour_ref_count := NA]
      region_bam_read_counts[line_idx, tumour_alt_count := NA]
      region_bam_read_counts[line_idx, tumour_N_count := NA]
    }
  }

  tumour_bam_read_counts <- rbindlist(list(tumour_bam_read_counts, region_bam_read_counts))

}

# No germline BAMs to count reads from -- every germline_*_count field is NA by construction.
tumour_bam_read_counts[, germline_ref_count := NA_real_]
tumour_bam_read_counts[, germline_alt_count := NA_real_]
tumour_bam_read_counts[, germline_N_count := NA_real_]
tumour_bam_read_counts[, germline_sample_name := NA_character_]

# add back the VEP CSQs etc (identical logic to upstream)
tumour_mutations[, c("mut_type", "start_dna", "end_dna") := NULL]
tumour_bam_read_counts[, c("start_dna", "end_dna") := NULL]
tumour_bam_read_counts <- merge(tumour_bam_read_counts, tumour_mutations,
                                by = c("CHROM", "POS", "REF", "ALT", "tumour_sample_name"),
                                all.x = TRUE)

# add in if it was called by mutect (identical logic to upstream)
tumour_bam_read_counts[is.na(FILTER), called_by_mutect := FALSE]
tumour_bam_read_counts[!is.na(FILTER), called_by_mutect := TRUE]

setnames(tumour_bam_read_counts,
         c("CHROM", "POS", "REF", "ALT", "FILTER", "CSQ"),
         c("chrom", "pos", "ref", "alt", "mutect_filter", "vep_csq"))

fwrite(tumour_bam_read_counts, file = mutation_save_path)
