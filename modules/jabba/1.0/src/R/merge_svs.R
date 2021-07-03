suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(purrr)
  library(gGnome)
  #library(StructuralVariantAnnotation)
  #library(GenomeInfoDb)
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
# args$manta <- '/projects/nhl_meta_analysis_scratch/gambl/results_local/gambl/manta-2.3/03-augment_vcf/genome--hg38/01-15563T--01-15563N--matched/somaticSV.augmented.vcf'
# args$gridss <- '/projects/nhl_meta_analysis_scratch/gambl/results_local/gambl/gridss-1.0/04-gripss/genome--hg38/01-15563T--01-15563N--matched/gridss_somatic_filtered.vcf.gz'
# args$genome <- 'hg38'
# args$pad <- 50

##### FUNCTIONS #####

# removeNonCanonicalChrs <- function(sv.vcf, chr_prefixed = TRUE, include_y = FALSE) {
#   chrs <- c(1:22, 'X')
#   if (include_y) {
#     chrs <- c(chrs, 'Y')
#   }
#   if (chr_prefixed) {
#     chrs <- paste0('chr', chrs)
#   }
#   
#   noncan <- sv.vcf[!seqnames(sv.vcf) %in% chrs]
#   removeIDs <- elementMetadata(noncan)
#   removeIDs <- unique(c(removeIDs[, 'sourceId'], removeIDs[, 'partner']))
#   
#   sv.vcf.filt <- sv.vcf[!names(sv.vcf) %in% removeIDs]
#   
#   return(sv.vcf.filt)
# }

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

gr2grl <- function(gr, genome) {
  nm  <- names(gr)
  out <- list()
  si <- Seqinfo(genome = genome)
  
  while (length(nm)) {
    s <- nm[1]
    p <- gr[s]$partner
    nm <- nm[!nm %in% c(s,p)]
    #out[[length(out)+1]] <- sort(c(gr.start(gr[s]), gr.start(gr[p])))
    out[[length(out)+1]] <- sort(c(gr[s], gr[p]))
  }
  
  return(GRangesList(out))
}

dt2bedpe <- function(dt) {
  dt %>% split(by = 'grl.ix') %>%
    map(function(x) {
      bp1 <- x[1,]; bp2 <- x[2,]
      data.table(chrom1 = bp1$seqnames, start1 = bp1$start, end1 = bp1$end, 
                 chrom2 = bp2$seqnames, start2 = bp2$start, end2 = bp2$end,
                 name = bp1$grl.ix, 
                 score = bp1$tier, 
                 strand1 = bp1$strand, 
                 strand2 = bp2$strand)
    }) %>%
    rbindlist()
}

##### MERGING #####

br <- list()
br$manta <- jJ(args$manta)[FILTER == 'PASS']
br$gridss <- jJ(args$gridss)[FILTER == 'PASS']

br <- map(br, ~ removeNonCanonicalChrs(.x))
br.merged <- merge(manta = br$manta, gridss = br$gridss, pad = args$pad)

br.merged.grl <- br.merged$grl
br.merged.val <- values(br.merged.grl)
br.merged.val$tier <- br.merged.val$seen.by.manta * br.merged.val$seen.by.gridss
br.merged.val$tier <- ifelse(br.merged.val$tier == 0, 2, 1)
values(br.merged.grl) <- br.merged.val

# br <- list()
# br$manta  <- breakpointRanges(readVcf(args$manta, genome = args$genome))
# br$gridss <- breakpointRanges(readVcf(args$gridss, genome = args$genome))
# 
# br <- map(br, function(x) x[x$FILTER == 'PASS',])
# br <- map(br, ~ removeNonCanonicalChrs(.x))
# 
# br.grl <- map(br, ~ gr2grl(.x, args$genome))
# br.merged <- grl.unlist(ra.merge(manta = br.grl$manta, gridss = br.grl$gridss, pad = args$pad))

# br.meta <- as.data.table(elementMetadata(br.merged))[, tier := 2][(seen.by.manta) & (seen.by.gridss), tier := 1]
# elementMetadata(br.merged) <- br.meta

# br.merged.dt <- gr2dt(br.merged)[,!c("grl.iix")]
# br.merged.dt[, tier := 2]
# br.merged.dt[(seen.by.manta) & (seen.by.gridss), tier := 1]

##### OUTPUT #####

saveRDS(br.merged.grl, file = args$rds)

# output <- GRangesList(split(gr2dt(br.merged), by = 'grl.ix'))
# saveRDS(output, file = args$rds)

# bedpe.out <- dt2bedpe(br.merged.dt)
# fwrite(bedpe.out, file = "/dev/stdout", sep = '\t', quote = FALSE, col.names = FALSE)








