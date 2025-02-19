#!/usr/bin/env Rscript

# Script to run MutationTimeR on one sample, inteded to be used in /projects/rmorin/projects/tumor-timing/src/MutationTimeR.smk
# Right now it assumes the sample has been run through Battenberg. It uses the filled subclones.txt (native or lifted) and
#   the cellularity_ploidy.txt for the purity estimate
# It assumes the tumour sample's data is in the merged, augmented, deblacklisted maf and uses the tumour_id to subset it on input


# Usage:
#   Rscript MutationTimeR.R <path/to/battenberg/subclones.txt> <path/to/battenberg/cellularity.txt> <path/to/augmented.maf> <path/to/output/SSM.tsv> <path/to/output/CNA.tsv> <tumour_sample_id> <projection> <path/for/log>

# Load packages -----------------------------------------------------------
cat("Loading packages...\n")
suppressWarnings(
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(MutationTimeR)
})
)

select = dplyr::select
filter = dplyr::filter
count = dplyr::count

# Parse command-line arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("bb_path", "cellularity_path", "maf_path", "n_bootstrap", "output_ssm_path", "output_cna_path", "tumour_sample_id", "projection", "log_path")
args <- setNames(args, arg_names[1:length(args)])

# Log both the stdout and stderr
log <- file(args$log_path, open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Print args for de-bugging -----------------------------------------------------------
cat(paste("bb_path:", args$bb_path, "\n"))
cat(paste("cellularity_path:", args$cellularity_path, "\n"))
cat(paste("augmented_maf_path:", args$maf_path, "\n"))
cat(paste("n_bootstrap:", args$n_bootstrap, "\n"))
cat(paste("output_ssm_path:", args$output_ssm_path, "\n"))
cat(paste("output_cna_path:", args$output_cna_path, "\n"))
cat(paste("tumour_sample_id:", args$tumour_sample_id, "\n"))
cat(paste("projection:", args$projection, "\n"))
cat(paste("log_path:", args$log_path, "\n"))

# Set command-line argument values -----------------------------------------------------
bb_file <- args$bb_path
cellularity_file <- args$cellularity_path
maf_file <- args$maf_path

n_bootstrap <- as.numeric(args$n_bootstrap)

output_ssm <- args$output_ssm_path
output_cna <- args$output_cna_path

tumour_id <- args$tumour_sample_id

this_projection <- args$projection

# Read in augmented maf data and convert to VCF-class
# -----------------------------------------------------
cat("Getting augmented maf SSM data...\n")
relevant_maf_columns <- c("Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome", "Start_Position", "End_Position", "Strand",
                          "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                          "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
                          "Transcript_ID", "Exon_Number", "t_depth", "t_ref_count", "t_alt_count","n_depth", "n_ref_count", "n_alt_count")

grep_cmd <- paste0('egrep "Tumor_Sample_Barcode|', tumour_id,'" ', maf_file, '| cut -f1-45')
maf <- fread(cmd = grep_cmd, verbose = F, nThread=4) %>%
  as_tibble() %>%
  select(all_of(relevant_maf_columns)) %>%
  mutate(Chromosome = as.character(Chromosome))

if(dim(maf)[1] == 0) stop(paste("Tumour sample", tumour_id, "is not in the input maf\n"))

cat("Subsetting to those with min read support at least 3...\n")
maf <- maf %>%
  filter(t_alt_count >= 3)

# TO DO: remove this filter when the vcf2maf deblacklist step has been fixed
maf <- unique(maf)

# chr prefix in the SSM data will be used to modify the CNA input so that they match
prefix_status <- ifelse(str_detect(maf$Chromosome[1], "chr"), TRUE, FALSE)

# Turn SSMs into VCF-class object
# -----------------------------------------------------
tumour_rowRanges <- maf %>%
  select(Chromosome, Start_Position, End_Position, Strand, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2) %>%
  makeGRangesFromDataFrame(., keep.extra.columns=TRUE, seqnames.field="Chromosome", start.field="Start_Position", end.field="End_Position",
    strand.field="Strand")
# make the DataFrame for the info header
info_descriptions <- data.frame(
  Number = c(1,1,1,1,1,1),
  Type = c("Integer","Integer","Integer","Integer","Integer","Integer"),
  Description = c("Depth in tumour sample", "Tumour ref count", "Tumour alt count", "Depth in normal sample", "Normal ref count", "Normal alt count"),
  row.names = c("t_depth", "t_ref_count", "t_alt_count", "n_depth", "n_ref_count", "n_alt_count")) %>%
  DataFrame(.)
# adding it to a blank VCFHeader object
header <- VCFHeader()
info(header) <- info_descriptions
# this makes the info DataFrame object
tumour_info <- maf %>%
  select(t_depth, t_ref_count, t_alt_count, n_depth, n_ref_count, n_alt_count) %>%
  DataFrame(.)
# this makes the VCF-class object
tumour_vcf <- VCF(rowRanges = tumour_rowRanges, exptData = list(header = header), info = tumour_info, collapsed = FALSE)


# Read in and format Battenberg output -----------------------------------------------------
cat("Reading in Battenberg files and formatting...\n")
purity <- read_tsv(cellularity_file, show_col_types=FALSE) %>% pull(cellularity)

bb <- read_tsv(bb_file, show_col_types=FALSE, na="NA") %>%
  dplyr::select(chr, startpos, endpos, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A)

# Filtering regions with nMaj1_A NA and nMin1_A NA -----------------------------------------------------------
# These are due to an oddity in the battenberg results, root cause has not been found yet
# Fixing chromosome prefixes to match maf -----------------------------------------------------------
# Fix the filled segments subclonal states -----------------------------------------------------------
# This is temporary until the battenberg results are regenerated
bb <- bb %>%
  filter(!(is.na(nMaj1_A) & is.na(nMin1_A))) %>%
  mutate(chr = case_when(
    prefix_status==TRUE & str_detect(chr, "chr") ~ chr, # both prefixed already, no changes needed
    prefix_status==FALSE & !str_detect(chr, "chr") ~ chr, # both not prefixed, no changes needed
    prefix_status==TRUE & !str_detect(chr, "chr") ~ paste0("chr", chr), # SSM prefixed, bb not, needs to be added
    prefix_status==FALSE & str_detect(chr, "chr") ~ gsub("chr", "", chr), # SSM not prefixed, bb is, needs to be removed
  )) %>%
  mutate(
    nMaj2_A = ifelse(frac1_A == 1 & frac2_A == 1, NA, nMaj2_A),
    nMin2_A = ifelse(frac1_A == 1 & frac2_A == 1, NA, nMin2_A),
    frac2_A = ifelse(frac1_A == 1 & frac2_A == 1, NA, frac2_A)
  )

bb_clonal <- bb %>%
  mutate(clonal_frequency = frac1_A*purity) %>%
  dplyr::select(chr, startpos, endpos, nMaj1_A, nMin1_A, clonal_frequency) %>%
  dplyr::rename(seqnames=chr, start=startpos, end=endpos, major_cn=nMaj1_A, minor_cn=nMin1_A)

bb_subclonal <- bb %>%
  filter(!is.na(frac2_A)) %>%
  mutate(clonal_frequency = frac2_A*purity) %>%
  dplyr::select(chr, startpos, endpos, nMaj2_A, nMin2_A, clonal_frequency) %>%
  dplyr::rename(seqnames=chr, start=startpos, end=endpos, major_cn=nMaj2_A, minor_cn=nMin2_A)

bb_final <- bind_rows(bb_clonal, bb_subclonal) %>% arrange(seqnames, start) %>% GRanges()

# Run MutationTimeR -----------------------------------------------------
cat("Running MutationTimeR...\n")
mt <-  mutationTime(tumour_vcf, bb_final, n.boot=n_bootstrap)

tumour_vcf <- addMutTime(tumour_vcf, mt$V)
mcols(bb_final) <- cbind(mcols(bb_final), mt$T)

# Write out results -----------------------------------------------------
# Check if output dir extists, create if not
output_dir <- dirname(output_cna)
if(!dir.exists(file.path(output_dir))){
  dir.create(file.path(output_dir), recursive=TRUE)
}

cat("Writing timed SSM outputs...\n")
# Get the info section as tibble
mt_vcf_tibble <- info(tumour_vcf) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(CNID_unnest = toString(CNID) %>% na_if(""), pAllSubclones_unnest = toString(pAllSubclones) %>% na_if("")) %>%
  ungroup() %>%
  select(-CNID, -pAllSubclones) %>%
  dplyr::rename(CNID=CNID_unnest, pAllSubclones=pAllSubclones_unnest)
# Get the rowRanges section
mt_vcf_rowRanges <- rowRanges(tumour_vcf)
# Add info section tibble to rowRanges
mcols(mt_vcf_rowRanges) <- cbind(mcols(mt_vcf_rowRanges), mt_vcf_tibble)
# Get final tibble of rowRanges data and info data
# Change CLS colname to more descriptive so users know it has the time labels
mt_vcf_tibble_final <- as_tibble(mt_vcf_rowRanges) %>%
  dplyr::rename(Chromosome=seqnames, Start_Position=start, End_Position=end, Strand=strand, CLS_time_label=CLS) %>%
  select(-width, -REF, -ALT, -QUAL, -FILTER) %>%
  relocate(CLS_time_label, .after = last_col())

# Join back to GAMBLR SSM data
mt_vcf_tibble_final_full <- inner_join(maf, mt_vcf_tibble_final, relationship = "one-to-one")

write_tsv(mt_vcf_tibble_final_full, file=output_ssm)

cat("Writing timed CNA outputs...\n")

bb_final_tibble <- as_tibble(bb_final) %>%
  select(-timing_param, -strand) %>%
  dplyr::rename(chr=seqnames, startpos=start, endpos=end)

write_tsv(bb_final_tibble, file=output_cna)

cat("DONE!")
sink()