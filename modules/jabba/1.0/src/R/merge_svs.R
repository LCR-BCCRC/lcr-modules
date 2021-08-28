suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(purrr)
  library(gGnome)
})

#####   PARSER  #####

parser <- ArgumentParser()
parser$add_argument('--manta', help = 'Path to Manta VCF.')
parser$add_argument('--gridss', help = 'Path to GRIDSS VCF.')
parser$add_argument('--genome', help = 'Genome build.')
parser$add_argument('--pad', default = 100, help = 'Padding for SV merging. [Default %(default)s].')
parser$add_argument('--rds', help = 'Output RDS')
args <- parser$parse_args()

##### TESTING #####

# args <- list()
# args$manta <- '/projects/nhl_meta_analysis_scratch/gambl/results_local/gambl/manta-2.3/03-augment_vcf/genome--grch37/16-20912T--16-20912N--matched/somaticSV.augmented.vcf'
# args$gridss <- '/projects/nhl_meta_analysis_scratch/gambl/results_local/gambl/gridss-1.0/04-gripss/genome--grch37/16-20912T--16-20912N--matched/gridss_somatic_filtered.vcf.gz'
# args$genome <- 'grch37'
# args$pad <- 50

##### FUNCTIONS #####

removeNonCanonicalChrs <- function(junc, chr_prefixed = TRUE, include_y = FALSE) {
  chrs <- c(1:22, 'X')
  if (include_y) {
    chrs <- c(chrs, 'Y')
  }
  if (chr_prefixed) {
    chrs <- paste0('chr', chrs)
  }
  
  junc <- junc[CHROM %in% chrs & MATECHROM %in% chrs]

  return(junc)
}

##### MERGING #####

br <- list()

sl <- NULL
if (args$genome == 'grch37') {
    sl <- hg_seqlengths('BSgenome.Hsapiens.UCSC.hg19::Hsapiens', chr = FALSE)
} else if (args$genome == 'hg38') {
    sl <- hg_seqlengths('BSgenome.Hsapiens.UCSC.hg38::Hsapiens', chr = TRUE)
}

br$manta <- jJ(args$manta, seqlengths = sl)
br$gridss <- jJ(args$gridss, seqlengths = sl)

#if (length(br$manta) != 0) {
#    br$manta <- br$manta[FILTER == 'PASS']
#}
#if (length(br$gridss) != 0) {
#    br$gridss <- br$gridss[FILTER == 'PASS']
#}

#if (any(args$genome %in% c("grch38", "grch37"))) {
#    br <- map(br, ~ removeNonCanonicalChrs(.x, chr_prefixed = FALSE))
#} else {
#    br <- map(br, ~ removeNonCanonicalChrs(.x))
#}

#empty.junc <- names(br)[which(lengths(br) == 0)]
br.merged <- merge(manta = br$manta, gridss = br$gridss, pad = args$pad)
br.merged.val <- as.data.table(values(br.merged$grl))
br.merged.grl <- GRangesList(br.merged$grl)
#br.merged.grl <- GRangesList(mapply(trim, br.merged$grl))
#br.merged.grl <- GRangesList(mapply(function(x) {ranges(x[start(x)<0] <- IRanges(1))}, br.merged$grl))

#for (i in empty.junc) {
#    br.merged.val[[i]] <- as.integer(rep(NA, nrow(br.merged.val)))
#    br.merged.val[[paste0('seen.by.',i)]] <- rep(FALSE, nrow(br.merged.val))
#}

# Categorize junctions into tiers
br.merged.val$tier <- br.merged.val$seen.by.manta * br.merged.val$seen.by.gridss
br.merged.val[, tier := ifelse(tier == 0, 2, 1)]
br.merged.val[FILTER.manta != "PASS" | FILTER.gridss != "PASS", tier := 3]
values(br.merged.grl) <- br.merged.val

##### OUTPUT #####

saveRDS(br.merged.grl, file = args$rds)

