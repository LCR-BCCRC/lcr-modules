#!/usr/bin/env Rscript
#
# Usage:
#   CNVfilteR_prepare.R <input_seg_file> <output_file> <sample_id>
#


suppressWarnings(suppressPackageStartupMessages({
  
  message("Loading packages...")
  
  required_packages <- c(
    "BiocManager", "dplyr", "data.table", "gtools", "readr", "tidyr", "tidyverse"
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
    "CNVfilteR", "rtracklayer", "BSgenome.Hsapiens.UCSC.hg38.masked", "InteractionSet"
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

SAMPLE.ID <- snakemake@wildcards[["tumour_sample_id"]]

SAMPLE.SEG <- filled.SEG %>% dplyr::filter(sample==SAMPLE.ID) %>%
                              dplyr::filter(!between(log.ratio,-0.4,0.4)) %>%
                              dplyr::mutate(sample="TUMOR") %>%
                              dplyr::mutate(CNV.type=ifelse(log.ratio<0, "deletion", "duplication"))
colnames(SAMPLE.SEG) <- c("sample", "chrom", "start", "end", "N.BAF", "log.ratio", "CNV.type")
# Output data for CNVfilteR ----------------------------------------------
message("Writing final table to file")
write_tsv(SAMPLE.SEG, snakemake@output[[1]])
Sys.sleep(5)

# Use CNVfilteR to filter CNVs --------------------------------------------
message("Loading data for CNV filtering...")

cnvs.file <- as.character(snakemake@output[[1]])

cnvs.gr <- loadCNVcalls(cnvs.file = cnvs.file, chr.column = "chrom", start.column = "start", end.column = "end", cnv.column = "CNV.type", sample.column = "sample", genome = "hg38")
cnvs.gr <- trim(cnvs.gr)
blacklist <-  import(snakemake@input[[3]])

vcf <- snakemake@output[[2]]
vcfs <- loadVCFs(vcf, cnvs.gr = cnvs.gr, vcf.source = "SLMS-3", list.support.field = "AD", genome = "hg38", regions.to.exclude = blacklist, min.total.depth = 4)

message("Filtering CNVs...")
results <- filterCNVs(cnvs.gr, vcfs, dup.threshold.score = 0.2, ht.deletions.threshold = 1, margin.pct = 25)
filtered <- results$cnvs[results$cnvs$filter != TRUE]
filtered <- cbind(as.data.frame(c(data.frame(filtered@seqnames),data.frame(filtered@ranges))),as.data.frame(filtered$cnv))

message("Preparing output...")
before <- SAMPLE.SEG %>% dplyr::select(sample, chrom, start, end, N.BAF, log.ratio)
before$sample <- snakemake@wildcards[["tumour_sample_id"]]

output <- merge(filtered, before, by="end") %>% dplyr::select(sample, chrom, start.y, end, N.BAF, log.ratio)
colnames(output) <- c("sample", "chr", "start", "end", "N.BAF", "log.ratio")
colnames(before) <- colnames(output)

# Keep large CNVs even if they are not supported by SNVs ------------
filtered.CNV <- dplyr::setdiff(before, output)
MIN.WIDTH <- 5000000
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
  mutate(N.BAF=ifelse(is.na(N.BAF), toString(2), N.BAF)) %>%
  mutate(log.ratio=ifelse(is.na(log.ratio), 0.0, log.ratio))

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

# sort output based on chromosome and segment start position
output <- output[multi.mixedorder(output$chr, output$start),]


# Output filtered data --------------------------------------------------
message("Writing final table to file")
write_tsv(output, snakemake@output[[2]])
