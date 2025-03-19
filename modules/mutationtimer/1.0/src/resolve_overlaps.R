#!/usr/bin/env Rscript

# Script to resolve overlapping segments in a Copy Number bed file made
#  from Battenberg subclones_filled.txt output that has been run through lcr-scripts/liftover.sh
# If the overlapping regions are clonal, 1:1 (without CN changes), they are merged
# If one of the overlapping regions has CNA info, that info is kept for the overlapping part
# If both segments have CNA info, the larger segment's info gets assigned to the overlapping region
#  the info from the smaller segment is written to another bed file, for the overlapping region
#  This bed file is named the same as the output_bed but with "_removed_in_ties.bed"
# Filters out regions in non-canonical chrs and where nMaj and nMin are clonally NA

#  Usage:
# Rscript <input.bed> <output.bed> <log.txt>

# Interactive:
# log_file <- "/projects/rmorin/projects/tumor-timing/test_withsameCNA.txt" # just as a placeholder, don't run the log related code interactively
# input_bed <- "/projects/rmorin/projects/tumor-timing/overlapping_withsameCNAs.bed"
# output_bed <- "/projects/rmorin/projects/tumor-timing/script_resolved_withsameCNAs.bed"
# args <- list(input_bed, output_bed, log_file)
# arg_names <- c("input_bed", "output_bed", "log_file")
# args <- setNames(args, arg_names[1:length(args)])



# Load packages -----------------------------------------------------------
suppressWarnings(
suppressPackageStartupMessages({
    library(tidyverse)
})
)

select = dplyr::select
filter = dplyr::filter
count = dplyr::count

# Parse command-line arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("input_bed", "output_bed", "log_file")
args <- setNames(args, arg_names[1:length(args)])

# Log both the stdout and stderr
log <- file(args$log_file, open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Print args for de-bugging -----------------------------------------------------------
removed_bed <- gsub("\\.bed$", "_removed_in_ties.bed", args$output_bed)
cat(paste("Input bed:", args$input_bed, "\n"))
cat(paste("Output bed:", args$output_bed, "\n"))
cat(paste("Bed file to save info from smaller segment in the case where both segments have CNA info:", removed_bed, "\n"))
cat(paste("log_file:", args$log_file, "\n"))


# Read in input and sort -----------------------------------------------------------
# Since it comes from lcr-scripts/liftover.sh, it should have a header

cat("Reading in and formatting input bed...\n")
bb_bed <- read_tsv(args$input_bed, show_col_types=FALSE, na="NA")

# Remove non-canonical chroms and chrY -----------------------------------------------------------
# Battenberg only reports on chrX, not chrY, so any chrY rows are due to the fill and are not meaningful
cat("Filtering non-canonical chromosomes and chrY...\n")

if (str_detect(bb_bed$chr[1], "chr")){
    chr_order <- c(paste0("chr",1:22),"chrX")
} else {
    chr_order <- c(1:22,"X")
}

bb_bed$chr <- factor(bb_bed$chr, levels=chr_order)

bb_bed <- bb_bed %>%
    dplyr::rename(chrom=chr, start=startpos, end=endpos) %>%
    filter(!is.na(chrom)) %>% # filters non-canon and Y
    arrange(chrom, start, end)


# Filtering regions with nMaj1_A NA and nMin1_A NA -----------------------------------------------------------
# These are due to an oddity in the battenberg results, root cause has not been found yet
cat("Removing regions with clonal CN states of NA...\n")
num_clonal <- dim(bb_bed %>% filter(is.na(nMaj1_A) & is.na(nMin1_A)))[1]
cat(paste0("This many segments had clonal NAs and will be removed: ", num_clonal, "\n"))
bb_bed <- bb_bed %>%
  filter(!(is.na(nMaj1_A) & is.na(nMin1_A)))

# Functions for resolving overlaps -----------------------------------------------------------
check_overlap <- function(bed) {
    highest_end = 0
    overlap <- c()
    for (i in 1:nrow(bed)) {
        if (i>1 && bed$chrom[i] == bed$chrom[i-1]) {
            if (bed$start[i] > highest_end) {
                overlap[i] = "NOToverlap"
            }else{
                overlap[i] = "overlap"
            }
            if (bed$end[i] > highest_end) {
                highest_end = bed$end[i]
            }
        } else {
            highest_end = bed$end[i]
            overlap[i] = "NOToverlap"
        }
    }
    bed <- bed %>%
        mutate(overlap_status = overlap)
    return(bed)
}

solve_overlap <- function(bed, nonnormal_removed_in_ties) {
  num_overlap = which(bed$overlap_status == "overlap")
  num_pre_overlap_sorted = (unique(sort(c(num_overlap-1,num_overlap))))
  non_overlap = bed[-num_pre_overlap_sorted,]
  bed <- bed[num_pre_overlap_sorted,]
  for (i in 1:nrow(bed)) {
    # as rows get removed, i counter can be greater than # rows in current bed
    if(i>dim(bed)[1]){
      break
    }
    if (bed$overlap_status[i] == "overlap") {
      if (bed$end[i] < bed$end[i-1]) {  # because of the sort,i is entirely within i-1
        if (is.na(bed$frac2_A[i]) & bed$nMaj1_A[i]==1 & bed$nMin1_A[i]==1){ # i is normal, can be removed
          bed <- bed[-c(i), ]
        } else if (is.na(bed$frac2_A[i-1]) & bed$nMaj1_A[i-1]==1 & bed$nMin1_A[i-1]==1){ # i is non-normal,  i-1 normal, split up i-1
          new_row1 <- data.frame(chrom = bed$chrom[i],
                              start = bed$start[i-1],
                              end = bed$start[i]-1,
                              nMaj1_A = bed$nMaj1_A[i-1],
                              nMin1_A = bed$nMin1_A[i-1],
                              frac1_A = bed$frac1_A[i-1],
                              nMaj2_A = bed$nMaj2_A[i-1],
                              nMin2_A = bed$nMin2_A[i-1],
                              frac2_A = bed$frac2_A[i-1],
                              overlap_status = "NOToverlap")
          new_row2 <- data.frame(chrom = bed$chrom[i],
                              start = bed$end[i]+1,
                              end = bed$end[i-1],
                              nMaj1_A = bed$nMaj1_A[i-1],
                              nMin1_A = bed$nMin1_A[i-1],
                              frac1_A = bed$frac1_A[i-1],
                              nMaj2_A = bed$nMaj2_A[i-1],
                              nMin2_A = bed$nMin2_A[i-1],
                              frac2_A = bed$frac2_A[i-1],
                              overlap_status = "NOToverlap")
          bed$overlap_status[i] <- "NOToverlap" # keeping i so change it's status
          bed <- bed[-c(i-1),]
          bed <- bind_rows(bed, new_row1, new_row2)
        } else if (identical(bed$nMaj1_A[i-1], bed$nMaj1_A[i]) & identical(bed$nMin1_A[i-1], bed$nMin1_A[i]) & identical(bed$nMaj2_A[i-1], bed$nMaj2_A[i]) & identical(bed$nMin2_A[i-1], bed$nMin2_A[i])){ # both are not normal, but the same, drop i
          bed <- bed[-c(i), ]
        } else { # both are not normal, not the same, remove i and record it
          nonnormal_removed_in_ties <- rbind(nonnormal_removed_in_ties, bed[c(i),])
          bed <- bed[-c(i),]
        }
      } else if (bed$start[i] == bed$start[i-1]){ # same start, bc of the sort, i is the larger segment here
        if(is.na(bed$frac2_A[i-1]) & bed$nMaj1_A[i-1]==1 & bed$nMin1_A[i-1]==1){ # i-1 is normal, can be removed
          bed$overlap_status[i] <- "NOToverlap"# keeping i so change it's status
          bed <- bed[-c(i-1),]
        } else if (is.na(bed$frac2_A[i]) & bed$nMaj1_A[i]==1 & bed$nMin1_A[i]==1){ # i is normal but i-1 is not, split i
          new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$end[i-1]+1,
                                    end = bed$end[i],
                                    nMaj1_A = bed$nMaj1_A[i],
                                    nMin1_A = bed$nMin1_A[i],
                                    frac1_A = bed$frac1_A[i],
                                    nMaj2_A = bed$nMaj2_A[i],
                                    nMin2_A = bed$nMin2_A[i],
                                    frac2_A = bed$frac2_A[i],
                                    overlap_status = "NOToverlap")
          bed <- bed[-c(i),]
          bed <- bind_rows(bed, new_row1)
        } else if (identical(bed$nMaj1_A[i-1], bed$nMaj1_A[i]) & identical(bed$nMin1_A[i-1], bed$nMin1_A[i]) & identical(bed$nMaj2_A[i-1], bed$nMaj2_A[i]) & identical(bed$nMin2_A[i-1], bed$nMin2_A[i])){ # both non-normal, but the same, drop i-1
          bed$overlap_status[i] <- "NOToverlap"# keeping i so change it's status
          bed <- bed[-c(i-1),]
        } else { # both are non-normal but not the same, drop i-1 but record it
          nonnormal_removed_in_ties <- rbind(nonnormal_removed_in_ties, bed[c(i-1),])
          bed$overlap_status[i] <- "NOToverlap"# keeping i so change it's status
          bed <- bed[-c(i-1),]
        }
      } else { # same end so segment i is smaller, or both have unique parts
        if ((is.na(bed$frac2_A[i]) & bed$nMaj1_A[i]==1 & bed$nMin1_A[i]==1) & (is.na(bed$frac2_A[i-1]) & bed$nMaj1_A[i-1]==1 & bed$nMin1_A[i-1]==1)){ # both are normal, merge into one section
          new_row1 <- data.frame(chrom = bed$chrom[i],
                                      start = bed$start[i-1],
                                      end = bed$end[i],
                                      nMaj1_A = bed$nMaj1_A[i],
                                      nMin1_A = bed$nMin1_A[i],
                                      frac1_A = bed$frac1_A[i],
                                      nMaj2_A = bed$nMaj2_A[i],
                                      nMin2_A = bed$nMin2_A[i],
                                      frac2_A = bed$frac2_A[i],
                                      overlap_status = "NOToverlap")
          bed <- bed[-c(i-1, i),]
          bed <- bind_rows(bed, new_row1)
        } else if (is.na(bed$frac2_A[i-1]) & bed$nMaj1_A[i-1]==1 & bed$nMin1_A[i-1]==1){ # i must not be bc of first test
          # keep i and unique part of i-1
          new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$start[i-1],
                                    end = bed$start[i]-1,
                                    nMaj1_A = bed$nMaj1_A[i-1],
                                    nMin1_A = bed$nMin1_A[i-1],
                                    frac1_A = bed$frac1_A[i-1],
                                    nMaj2_A = bed$nMaj2_A[i-1],
                                    nMin2_A = bed$nMin2_A[i-1],
                                    frac2_A = bed$frac2_A[i-1],
                                    overlap_status = "NOToverlap")
          bed$overlap_status[i] <- "NOToverlap"
          bed <- bed[-c(i-1),]
          bed <- bind_rows(bed, new_row1)
        } else if(is.na(bed$frac2_A[i]) & bed$nMaj1_A[i]==1 & bed$nMin1_A[i]==1){ # i-1 must not be
          # keep i-1 and unqiue part of i, unless they had same ends
          new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$end[i-1]+1,
                                    end = bed$end[i],
                                    nMaj1_A = bed$nMaj1_A[i],
                                    nMin1_A = bed$nMin1_A[i],
                                    frac1_A = bed$frac1_A[i],
                                    nMaj2_A = bed$nMaj2_A[i],
                                    nMin2_A = bed$nMin2_A[i],
                                    frac2_A = bed$frac2_A[i],
                                    overlap_status = "NOToverlap")
          if(bed$end[i] == bed$end[i-1]){
            bed <- bed[-c(i),]
          } else {
            bed <- bed[-c(i),]
            bed <- bind_rows(bed, new_row1)
          }
        } else if (identical(bed$nMaj1_A[i-1], bed$nMaj1_A[i]) & identical(bed$nMin1_A[i-1], bed$nMin1_A[i]) & identical(bed$nMaj2_A[i-1], bed$nMaj2_A[i]) & identical(bed$nMin2_A[i-1], bed$nMin2_A[i])){ # both non-normal, but the same, merge into one
          new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$start[i-1],
                                    end = bed$end[i],
                                    nMaj1_A = bed$nMaj1_A[i],
                                    nMin1_A = bed$nMin1_A[i],
                                    frac1_A = bed$frac1_A[i],
                                    nMaj2_A = bed$nMaj2_A[i],
                                    nMin2_A = bed$nMin2_A[i],
                                    frac2_A = bed$frac2_A[i],
                                    overlap_status = "NOToverlap")
          bed <- bed[-c(i-1, i),]
          bed <- bind_rows(bed, new_row1)
        } else { # both have info but not the same
          if((bed$end[i]-bed$start[i]) > (bed$end[i-1]-bed$start[i-1])){
            # keep i and unique part of i-1, write out intersection's i-1 info
            new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$start[i-1],
                                    end = bed$start[i]-1,
                                    nMaj1_A = bed$nMaj1_A[i-1],
                                    nMin1_A = bed$nMin1_A[i-1],
                                    frac1_A = bed$frac1_A[i-1],
                                    nMaj2_A = bed$nMaj2_A[i-1],
                                    nMin2_A = bed$nMin2_A[i-1],
                                    frac2_A = bed$frac2_A[i-1],
                                    overlap_status = "NOToverlap")
            intersection <- data.frame(chrom = bed$chrom[i],
                                      start = bed$start[i],
                                      end = bed$end[i-1],
                                      nMaj1_A = bed$nMaj1_A[i-1],
                                      nMin1_A = bed$nMin1_A[i-1],
                                      frac1_A = bed$frac1_A[i-1],
                                      nMaj2_A = bed$nMaj2_A[i-1],
                                      nMin2_A = bed$nMin2_A[i-1],
                                      frac2_A = bed$frac2_A[i-1],
                                      overlap_status = "NOToverlap")
            nonnormal_removed_in_ties <- rbind(nonnormal_removed_in_ties, intersection)
            bed$overlap_status[i] <- "NOToverlap"
            bed <- bed[-c(i-1),]
            bed <- bind_rows(bed, new_row1)
          } else {
            # keep i-1 and unique part of i (unless they had the same end), write out intersection's i info
            new_row1 <- data.frame(chrom = bed$chrom[i],
                                    start = bed$end[i-1]+1,
                                    end = bed$end[i],
                                    nMaj1_A = bed$nMaj1_A[i],
                                    nMin1_A = bed$nMin1_A[i],
                                    frac1_A = bed$frac1_A[i],
                                    nMaj2_A = bed$nMaj2_A[i],
                                    nMin2_A = bed$nMin2_A[i],
                                    frac2_A = bed$frac2_A[i],
                                    overlap_status = "NOToverlap")
            intersection <- data.frame(chrom = bed$chrom[i],
                                      start = bed$start[i],
                                      end = bed$end[i-1],
                                      nMaj1_A = bed$nMaj1_A[i],
                                      nMin1_A = bed$nMin1_A[i],
                                      frac1_A = bed$frac1_A[i],
                                      nMaj2_A = bed$nMaj2_A[i],
                                      nMin2_A = bed$nMin2_A[i],
                                      frac2_A = bed$frac2_A[i],
                                      overlap_status = "NOToverlap")
            nonnormal_removed_in_ties <- rbind(nonnormal_removed_in_ties, intersection)
            if(bed$end[i] == bed$end[i-1]){
              bed <- bed[-c(i),]
            } else {
              bed <- bed[-c(i),]
              bed <- bind_rows(bed, new_row1)
            }
          }
        }
      }
    }
    bed <- bed %>%
      arrange(chrom,start,end)
  }

  bed <- bind_rows(bed, non_overlap) %>%
    arrange(chrom,start,end)

  nonnormal_removed_in_ties <- nonnormal_removed_in_ties %>%
    arrange(chrom,start,end)

  bed <- check_overlap(bed)
  while("overlap" %in% bed$overlap_status){
    solve_list <- solve_overlap(bed, nonnormal_removed_in_ties)
    bed <- solve_list[[1]]
    nonnormal_removed_in_ties <- solve_list[[2]]
  }

  return(list(bed, nonnormal_removed_in_ties))
}

# Actual resolving overlaps -----------------------------------------------------------
cat("Checking for overlaps...\n")
bb_bed_checked <- check_overlap(bb_bed)

# df to save the removed non-normal segments that lost ties
cols <- colnames(bb_bed_checked)
removed_in_ties <- data.frame(matrix(nrow=0, ncol=length(cols)))
colnames(removed_in_ties) <- cols

# Only send through solve function if there are overlaps
if (sum(bb_bed_checked$overlap_status == "overlap") == 0){
  cat("No overlaps detected, writing outputs...\n")
  # Check if output dir extists, create if not
  output_dir <- dirname(args$output_bed)
  if(!dir.exists(file.path(output_dir))){
    dir.create(file.path(output_dir), recursive=TRUE)
  }

  bb_bed_resolved <- bb_bed %>%
      dplyr::rename(chr=chrom, startpos=start, endpos=end)

  write_tsv(bb_bed_resolved, file=args$output_bed)

  write_tsv(removed_in_ties, file=removed_bed)
} else {
  cat("Resolving overlaps...\n")
  solve_overlaps_list <- solve_overlap(bb_bed_checked, removed_in_ties)
  bb_bed_resolved <- solve_overlaps_list[[1]]
  removed_in_ties <- solve_overlaps_list[[2]]

  cat("Writing outputs...\n")
  # Check if output dir extists, create if not
  output_dir <- dirname(args$output_bed)
  if(!dir.exists(file.path(output_dir))){
    dir.create(file.path(output_dir), recursive=TRUE)
  }

  bb_bed_resolved <- bb_bed_resolved %>%
      dplyr::rename(chr=chrom, startpos=start, endpos=end)

  write_tsv(bb_bed_resolved, file=args$output_bed)

  write_tsv(removed_in_ties %>% select(-overlap_status), file=removed_bed)
}

cat("DONE!")
sink()