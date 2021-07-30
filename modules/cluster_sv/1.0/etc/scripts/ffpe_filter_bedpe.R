# /usr/bin/Rscript
# Author: Helena Winata 
# adapted from https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R
#################################################################################################### 
# USAGE: 
#     size_filter_BEDPE.R -bedpe <input.bedpe> -svtype <SVTYPE> -thres <size.thres> -out_bedpe <output.bedpe>
# Arguments:
#   input.bedpe (-bedpe) - Input BEPDE file with annotation or SV_ID in the format of INV_<id>.
#   SVTYPE (-svtype) - SVTYPE to filter out (INV, DUP, DEL, INS, ITX, etc.); Default = INV
#   size.thres (-thres) - threshold for size to filter out; Default = 0
#   output.bedpe (-out_bedpe) - Optional output filename (default: <input_filename>_filt.bedpe)


library(R.utils)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(rtracklayer)
library(stringr)


# Read command line arguments
ARGS <- commandArgs(
  trailingOnly = FALSE,
  asValues = TRUE,
  defaults = list(ref = 'hg38', thres = 0, svtype = 'INV')
);


file <- ARGS$bedpe;
if (!file.exists(file)) {
  stop(sprintf('Input file %s does not exist!', file))
};

if (is.null(ARGS$out_bedpe)) {
  out_bedpe <- sub('.bedpe$', '_filt.bedpe', file)
} else {
  out_bedpe <- ARGS$out_bedpe
};

bedpe <- read.table(file = file, 
                    header = FALSE,
                    sep = '\t', 
                    col.names = c(
                      'CHROM_A', 'START_A', 'END_A',	
                      'CHROM_B', 'START_B', 'END_B',
                      'SV_ID', 'QUAL', 'STR_A', 'STR_B')
                    )

bedpe$SIZE <- abs(((bedpe$START_B + bedpe$END_B) / 2) - ((bedpe$START_A + bedpe$END_A) / 2))
bedpe$SVTYPE <- gsub('_.*', '', bedpe$SV_ID)

sprintf('Removing BPs with SVTYPE = %s with size less than %.0f', ARGS$svtype, ARGS$thres)
filt_bedpe <- subset(bedpe, SIZE >= ARGS$thres & SVTYPE != ARGS$svtype)

write.table(filt_bedpe[, c('CHROM_A', 'START_A', 'END_A',	
                            'CHROM_B', 'START_B', 'END_B',
                            'SV_ID', 'QUAL', 'STR_A', 'STR_B')],
            file = out_bedpe, sep = '\t', 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
