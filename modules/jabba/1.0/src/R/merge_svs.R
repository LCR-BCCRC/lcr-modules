suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(purrr)
  library(gGnome)
  library(StructuralVariantAnnotation)
})

#####   PARSER  #####

parser <- ArgumentParser()
parser$add_argument('--manta', type = , help = 'Path to Manta VCF.')
parser$add_argument('--gridss', help = 'Path to GRIDSS VCF.')
parser$add_argument('--pad', default = 100, help = 'Padding for SV merging. [Default %(default)s].')
args <- parser$parse_args()

##### FUNCTIONS #####

removeNonCanonicalChrs <- function(sv.vcf, chr_prefixed = TRUE, include_y = FALSE) {
  chrs <- c(1:22, 'X')
  if (include_y) {
    chrs <- c(chrs, 'Y')
  }
  if (chr_prefixed) {
    chrs <- paste0('chr', chrs)
  }
  
  noncan <- sv.vcf[!seqnames(sv.vcf) %in% chrs]
  removeIDs <- elementMetadata(noncan)
  removeIDs <- unique(c(removeIDs[, 'sourceId'], removeIDs[, 'partner']))
  
  sv.vcf.filt <- sv.vcf[!names(sv.vcf) %in% removeIDs]
  
  return(sv.vcf.filt)
}

gr2grl <- function(gr) {
  nm  <- names(gr)
  out <- list()
  
  while (length(nm)) {
    s <- nm[1]
    p <- gr[s]$partner
    nm <- nm[!nm %in% c(s,p)]
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
                 score = NA, 
                 strand1 = bp1$strand, 
                 strand2 = bp2$strand,
                 tier = bp1$tier)
    }) %>%
    rbindlist()
}

##### MERGING #####

br <- list()
br$manta  <- breakpointRanges(readVcf(args$manta, genome = 'hg38'))
br$gridss <- breakpointRanges(readVcf(args$gridss, genome = 'hg38'))

br <- map(br, function(x) x[x$FILTER == 'PASS',])
br <- map(br, ~ removeNonCanonicalChrs(.x))

br.grl <- map(br, gr2grl)
br.merged <- grl.unlist(ra.merge(manta = br.grl$manta, gridss = br.grl$gridss, pad = args$pad))

##### TIER #####

br.merged.dt <- gr2dt(br.merged)[,!c("grl.iix")]
br.merged.dt[, tier := 2]
br.merged.dt[(seen.by.manta) & (seen.by.gridss), tier := 1]

##### OUTPUT #####

bedpe.out <- dt2bedpe(br.merged.dt)
fwrite(bedpe.out, file = "/dev/stdout", sep = '\t', quote = FALSE, col.names = FALSE)

