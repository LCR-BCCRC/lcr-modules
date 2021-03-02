#!/usr/bin/env Rscript
#
# Usage:
#   CNVfilteR_prepare.R <input_seg_file> <output_file> <sample_id>
#


suppressWarnings(suppressPackageStartupMessages({
  
  message("Loading packages...")
  
  required_packages <- c(
    "BiocManager", "dplyr", "data.table", "gtools", "readr", "tidyr", "tidyverse", "devtools"
  )
  
  installed_packages <- rownames(installed.packages())
  
  for (pkg in required_packages){
    if (!pkg %in% installed_packages) {
      message(paste0("Installing ", pkg, "..."))
      install.packages(pkg, repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
    }
    library(pkg, character.only = TRUE)
  }
  
  BiocManager_packages <- c(
    "CNVfilteR", "rtracklayer", "BSgenome.Hsapiens.UCSC.hg38.masked", "BSgenome.Hsapiens.UCSC.hg19.masked", "InteractionSet"
  )
  for (pkg in BiocManager_packages){
    if (!pkg %in% installed.packages()) {
      message(paste0("Installing ", pkg, " from GitHub..."))
      BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
  }  
}))


# Prepare data for CNVfilteR ----------------------------------------------
message("Loading data for pre-processing...")

filled.SEG <- fread(snakemake@input[[1]])
colnames(filled.SEG) <- c("sample", "chrom", "start", "end", "LOH_flag", "log.ratio")

SAMPLE.ID <- snakemake@wildcards[["tumour_id"]]

SAMPLE.SEG <- filled.SEG %>% dplyr::filter(sample==SAMPLE.ID) %>%
                              dplyr::filter(!between(log.ratio,-0.2,0.2)) %>%
                              dplyr::mutate(sample="TUMOR") %>%
                              dplyr::mutate(CNV.type=ifelse(log.ratio<0, "deletion", "duplication")) %>%
                              dplyr::mutate(chrom=ifelse(grepl("chr", chrom), chrom, paste0("chr",chrom)))

if(nrow(SAMPLE.SEG)<2) {
  ROW <- data.frame("sample" = "TUMOR",              
                  "chrom" = "chrX",
                  "start" = 1,
                  "end" = 300,
                  "LOH_flag" = 0,
                  "log.ratio" = 0.58,
                  CNV.type="duplication")
  SAMPLE.SEG <- rbind(SAMPLE.SEG, ROW)           
}


# Output data for CNVfilteR ----------------------------------------------
message("Writing final table to file ", snakemake@output[[1]])
write_tsv(SAMPLE.SEG, snakemake@output[[1]])
Sys.sleep(5)

# Use CNVfilteR to filter CNVs --------------------------------------------
message("Loading ", snakemake@output[[1]], " for CNV filtering...")

cnvs.file <- as.character(snakemake@output[[1]])
genome.build <- as.character(snakemake@params[["genome_build"]])
cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "chrom", start.column = "start", end.column = "end", cnv.column = "CNV.type", sample.column = "sample", genome = genome.build)
cnvs.gr <- trim(cnvs.gr)
blacklist <-  import(snakemake@input[[3]])

vcf <- snakemake@input[[2]]
vcfs <- loadVCFs(vcf, cnvs.gr = cnvs.gr, vcf.source = "SLMS-3", list.support.field = "AD", genome = genome.build, regions.to.exclude = blacklist, min.total.depth = 4)

message("Filtering CNVs...")
results <- filterCNVs(cnvs.gr, vcfs, dup.threshold.score = 0.2, ht.deletions.threshold = 1, margin.pct = 25)
filtered <- results$cnvs[results$cnvs$filter != TRUE]
filtered <- cbind(as.data.frame(c(data.frame(filtered@seqnames),data.frame(filtered@ranges))),as.data.frame(filtered$cnv))

message("Preparing output...")
before <- SAMPLE.SEG %>% dplyr::select(sample, chrom, start, end, LOH_flag, log.ratio)
before$sample <- snakemake@wildcards[["tumour_id"]]

output <- merge(filtered, before, by="end") %>% dplyr::select(sample, chrom, start.y, end, LOH_flag, log.ratio)
colnames(output) <- c("sample", "chr", "start", "end", "LOH_flag", "log.ratio")
colnames(before) <- colnames(output)

# Keep large CNVs even if they are not supported by SNVs ------------
filtered.CNV <- dplyr::setdiff(before, output)
MIN.WIDTH <- 4000000
if(nrow(filtered.CNV)>0) {
  for (row in 1:nrow(filtered.CNV)) {
    width <- as.numeric(filtered.CNV[row, "end"]) - as.numeric(filtered.CNV[row, "start"]) + 1
    
    if(width > MIN.WIDTH) {
      output <- rbind(output,filtered.CNV[row,])
    }
  }
}
chromosomes <- c()
for (x in 1:22) {
  chromosomes[x] = paste0("chr",x)
  if (x==22) {
    chromosomes[23] = "chrX"
  }
}

# make sure at least 1 entry is present for each chromosome
output <- complete(output, tidyr::expand(output, crossing(sample), chr = chromosomes))

output <- output %>% mutate(start=ifelse(is.na(start), paste0("100",rownames(.)), start)) %>%
  mutate(end=ifelse(is.na(end), paste0("100000",rownames(.)), end)) %>%
  mutate(LOH_flag=ifelse(is.na(LOH_flag), toString(2), LOH_flag)) %>%
  mutate(log.ratio=ifelse(is.na(log.ratio), 0.0, log.ratio)) %>%
  mutate(chr=ifelse(grepl("chr",chr), chr, paste0("chr",chr)))

# define function to use for mixed ordering based on multiple columns
multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

# add new segment att at the end to allow for segment filling
new <- data.frame("sample" = SAMPLE.ID,              
                  "chr" = "chrX",
                  "start" = 156040965,
                  "end" = 156040975,
                  "LOH_flag" = 2,
                  "log.ratio" = 0.0)
output <- rbind(output, new)

# sort output based on chromosome and segment start position
output <- output[multi.mixedorder(output$chr, output$start),]


# Output filtered data --------------------------------------------------
message("Writing final table to file")
write_tsv(output, snakemake@output[[2]])
