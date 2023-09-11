suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(purrr)
  library(gGnome)
  library(gUtils)
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

# setwd('/projects/rmorin/projects/gambl-repos/gambl-kcoyle/')
# args <- list()
# args$manta <- 'results/gambl/jabba_battenberg_purities/00-inputs/junc/manta/genome--hg38/01-20357T--01-20357N--matched.vcf'
# args$gridss <- 'results/gambl/jabba_battenberg_purities/00-inputs/junc/gridss/genome--hg38/01-20357T--01-20357N--matched.bnd.vcf'
# args$genome <- 'hg38'
# args$pad <- 100
# args$rds <- 'results/gambl/jabba_battenberg_purities/00-inputs/junc/merged/genome--hg38/01-20357T--01-20357N--matched.rds'

##### MERGING #####

br <- list()

sl <- NULL


if (args$genome == 'grch37') {
    sl <- hg_seqlengths('BSgenome.Hsapiens.UCSC.hg19::Hsapiens', chr = FALSE)
} else if (args$genome == 'hg38') {
    sl <- hg_seqlengths('BSgenome.Hsapiens.UCSC.hg38::Hsapiens', chr = TRUE)
}


# FIX: jJ fails when no breakpoints are present
# i.e., only del, dup, ins are present
# Should only be a problem with Manta

# Error arises when using seqlengths, removed ", seqlengths = sl" from manta & gridss
br$manta <- 
  tryCatch({
    jJ(args$manta#,  hg=args$genome)
    )
    #jJ(args$manta, seqlengths = hg_seqlengths(genome = DEFAULT_GENOME))
  }, error = function(e){
    GRanges()
  })
br$gridss <- jJ(args$gridss#,  hg=args$genome
               )

any.empty <- any(map_int(br, length)==0)

# Handle empty SV inputs
if (length(br$manta) == 0) {
  br.merged <- merge(gridss = br$gridss, pad = args$pad)
} else if (length(br$gridss) == 0) {
  br.merged <- merge(manta = br$manta, pad = args$pad)
} else if (length(br$gridss) == 0 && length(br$manta) == 0) {
  warning("GRIDSS and Manta VCFs are both empty.\n")
} else {
  br.merged <- merge(manta = br$manta, gridss = br$gridss, pad = args$pad)
}

br.merged.val <- as.data.table(values(br.merged$grl))
br.merged.grl <- GRangesList(br.merged$grl)

# Categorize variants into tiers if SVs are available from both tools
if (!any.empty) {
  br.merged.val$tier <- br.merged.val$seen.by.manta * br.merged.val$seen.by.gridss
  br.merged.val[, tier := ifelse(tier == 0, 2, 1)]
  br.merged.val[FILTER.manta != "PASS" | FILTER.gridss != "PASS", tier := 3]
}

values(br.merged.grl) <- br.merged.val

##### OUTPUT #####
if (!dir.exists(dirname(args$rds))) {dir.create(dirname(args$rds), recursive = TRUE)}
saveRDS(br.merged.grl, file = args$rds)

