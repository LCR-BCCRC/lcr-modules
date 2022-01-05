# /usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(stdout(), type = "message")
sink(log)

suppressWarnings({
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(StructuralVariantAnnotation)  
    library(data.table)
  })
  

  col_spec <- cols(
    CHROM_A = col_character(),
    START_A = col_double(),
    END_A = col_double(),
    CHROM_B = col_character(),
    START_B = col_double(),
    END_B = col_double(),
    manta_name = col_character(),
    SCORE = col_double(),
    STRAND_A = col_character(),
    STRAND_B = col_character(),
    tumour_sample_id = col_character(),
    normal_sample_id = col_character(),
    VAF = col_double(),
    DP = col_double(),
    gridss_name = col_character()
)
  
  ##### GET INPUT VALUES #####
  gridss <- snakemake@input[["gridss"]]
  manta <- snakemake@input[["manta"]]
  
  output_file <- snakemake@output[["bedpe"]]
  
  bed <- snakemake@params[["regions_bed"]]
  bedpe_cols <- snakemake@params[["bedpe_cols"]]
  minvaf <- snakemake@params[["minvaf"]]
  mindp <- snakemake@params[["mindp"]]
  maxgap <- snakemake@params[["maxgap"]]
  
  genome_build <- snakemake@wildcards[["genome_build"]]
  tumour_id <- snakemake@wildcards[["tumour_id"]]
  normal_id <- snakemake@wildcards[["normal_id"]]
  
  imagefile <- snakemake@log[[2]]
  
  cat(paste0("Starting SV filtering for tumour ", tumour_id, " on ", Sys.getenv("HOSTNAME")))
  
  ##### GET REGIONS BED FILE #####
  
  # Load the bed file describing the regions to whitelist
  cat("Loading the regions bed file...\n", quote = FALSE)
  
  
  bed <- read_tsv(bed,
                  col_names = c("chr", "start", "end", "name"),
                  col_types = cols()) %>%
    mutate(strand = "+")
  
  save.image(imagefile)
  
  bed_tbl <- data.table(bed)
  
  # Convert the bed file to GRanges Object
  
  bed <- with(bed, GRanges(chr, IRanges(start, end), strand, name))
  
  save.image(imagefile)
  
  # Define a function to identify all breakpoints with one or both ends in the capture space
  get_in_region <- function(vcf, bed){
    unique(sort(c(
      queryHits(findOverlaps(vcf, bed, ignore.strand = TRUE)), 
      queryHits(findOverlaps(partner(vcf, selfPartnerSingleBreakends = TRUE), bed, ignore.strand = TRUE)) 
    )))
  }
  
  # Define a function to whitelist all intrachromosomal variants between regions in the capture space
  get_in_whitelist <- function(vcf, bed){
    first_end <- vcf[queryHits(findOverlaps(vcf, bed, ignore.strand = TRUE))]
    second_end <- vcf[queryHits(findOverlaps(partner(vcf), bed, ignore.strand = TRUE))]
    first_end_matched <- first_end[first_end$partner %in% names(first_end)]
    intrachrom <- first_end_matched[seqnames(first_end_matched) != seqnames(partner(first_end_matched))]
    return(intrachrom)
  }
  
  
  # Define a function for depth and VAF filtering breakpoints
  filter_vaf_depth <- function(vcf, vaf_cut, depth_cut){
    vcf[(vcf$VAF >= vaf_cut & vcf$DP >= depth_cut) | (partner(vcf)$VAF >= vaf_cut & partner(vcf)$DP >= depth_cut)] 
  }
  
  save.image(imagefile)
  
  ##### GRIDSS #####
  
  # Parse the GRIDSS VCF file into a breakpoint ranges object
  
  cat("Loading GRIDSS VCF...\n")
  # gridss <-  breakpointRanges(readVcf("results/gridss-2.0/99-outputs/vcf/capture--hg38/CLC02071--capture_TE99028370--None--no_normal.gridss_somatic_filtered.vcf.gz",
  #                                     "hg38"),
  #                             info_columns = c("REF", "REFPAIR", "VF"))

  gridss <-  breakpointRanges(readVcf(gridss, 
                                      genome_build),
                              info_columns = c("REF", "REFPAIR", "VF"))
  
  # Calculate VAF and DP fields
  cat("Annotating GRIDSS VAF...")
  gridss_vaf <- elementMetadata(gridss)
  gridss_vaf$DP <- gridss_vaf$VF + gridss_vaf$REF + gridss_vaf$REFPAIR
  gridss_vaf$VAF <- gridss_vaf$VF/gridss_vaf$DP
  elementMetadata(gridss) <- gridss_vaf
  
  
  # Filter to variants with minimum VAF and depth
  cat("Filtering GRIDSS VCF...")
  gridss_filt <- filter_vaf_depth(gridss, minvaf, mindp)
  
  # Get whitelist breakpoints (not VAF filtered)
  gridss_whitelist <- get_in_whitelist(gridss, bed)
  gridss_whitelist <- gridss_whitelist[!names(gridss_whitelist) %in% names(gridss_filt)]
  
  # Combine whitelist with VAF-filtered breaks
  gridss_filt <- c(gridss_filt, gridss_whitelist)
  
  save.image(imagefile)
  
  ##### MANTA #####
  
  # Load the Manta vcf file
  cat("Loading Manta VCF...")
 

  manta_vcf <- readVcf(manta,
                       genome_build,
                       ScanVcfParam(geno = c("DP", "VAF")))
  
  # Put the Manta VAF and DP fields into the info field
  cat("Annotating Manta VAF...")
  info(manta_vcf)$VAF <- geno(manta_vcf)$VAF[ , tumour_id]
  info(manta_vcf)$DP <- geno(manta_vcf)$DP[ , tumour_id]
  
  # Convert the Manta VCF object to a breakpointRanges object
  manta <- breakpointRanges(manta_vcf, 
                            info_columns = c("DP", "VAF"))
  
  # Filter to only variants overlapping the capture space
  cat("Filtering Manta VCF...")
  manta_filt <- filter_vaf_depth(manta, minvaf, mindp)
  
  # Get whitelist breakpoints  (not VAF filtered)
  manta_whitelist <- get_in_whitelist(manta, bed)
  manta_whitelist <- manta_whitelist[!names(manta_whitelist) %in% names(manta_filt)]

  # Combine VAF filtered and whitelisted breakpoints
  manta_filt <- c(manta_filt, manta_whitelist)
  
  save.image(imagefile)
  
  ##### INTERSECT MANTA GRIDSS #####
  
  cat("Intersecting Manta and GRIDSS calls...\n")
  
  # Find the overlaps between Manta and GRIDSS
  overlaps_mg <- findBreakpointOverlaps(manta_filt, gridss_filt, maxgap = maxgap)
  
  overlaps_mg <- data.frame(
    gridss_name = names(gridss_filt[subjectHits(overlaps_mg)]), 
    manta_name = names(manta_filt[queryHits(overlaps_mg)])
  )
  
  save.image(imagefile)
  
  # Subeset the Manta vcf to those in the overlaps
  manta_overlaps <- manta_filt[unique(overlaps_mg[!is.na(overlaps_mg$manta_name), ]$manta_name)]
  
  # Convert to bedpe
  cat("Generating Manta BEDPE...\n")
  manta_bedpe <- breakpointgr2bedpe(manta_filt)
  # Update colnames and add some colums to match GAMBL standard
  colnames(manta_bedpe) <- c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SCORE",  "STRAND_A", "STRAND_B")
  if(nrow(manta_bedpe) > 0){
    manta_bedpe$tumour_sample_id <- tumour_id
    manta_bedpe$normal_sample_id <- normal_id
    manta_bedpe$VAF <- manta_filt[manta_bedpe$NAME, ]$VAF
    manta_bedpe$DP <- manta_filt[manta_bedpe$NAME, ]$DP
    manta_bedpe <- left_join(manta_bedpe, overlaps_mg, by = c("NAME" = "manta_name")) %>% 
      rename("NAME" = "manta_name")
  }
  
  save.image(imagefile)
  
  # Subset the GRIDSS vcf to those in the overlaps
  gridss_overlaps <- gridss_filt[unique(overlaps_mg[!is.na(overlaps_mg$gridss_name), ]$gridss_name)]
  
  # Convert to bedpe
  cat("Generating GRIDSS BEDPE...\n")
  gridss_bedpe <- breakpointgr2bedpe(gridss_filt)
  colnames(gridss_bedpe) <- c("CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "NAME", "SCORE",  "STRAND_A", "STRAND_B")
  if(nrow(gridss_bedpe) > 0){ 
    gridss_bedpe$tumour_sample_id <- tumour_id
    gridss_bedpe$normal_sample_id <- normal_id
    gridss_bedpe$VAF <- gridss_filt[gridss_bedpe$NAME, ]$VAF
    gridss_bedpe$DP <- gridss_filt[gridss_bedpe$NAME, ]$DP
    gridss_bedpe <- left_join(gridss_bedpe, overlaps_mg, by = c("NAME" = "gridss_name")) %>% 
      rename("NAME" = "gridss_name")
    
  }
  
  save.image(imagefile)
  
  # Combine Manta and GRIDSS bedpes and get unique rows based on variant name columns
  bedpe_combined <- read_tsv(bedpe_cols, 
                            col_types = col_spec
                    )
  bedpe_combined <- bind_rows(bedpe_combined, manta_bedpe, gridss_bedpe)
  
  if(nrow(overlaps_mg) > 0){
    
    bedpe_combined <- bedpe_combined  %>% 
    distinct(manta_name, gridss_name, .keep_all = TRUE)

    # Combine the overlaps bedpe with the Manta RNA bedpe if it exists

    
  
  }
  # Write combined bedpe to file
  write_tsv(dplyr::rename(bedpe_combined, SOMATIC_SCORE = SCORE), output_file)
  
  
  
  cat(paste0("SV filtering complete for tumour ", tumour_id, " on ", Sys.getenv("HOSTNAME")))
  
  save.image(imagefile)
  
})

sink()



