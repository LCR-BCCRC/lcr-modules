# /usr/bin/Rscript
# Author: Helena Winata 
# adapted from https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R
#################################################################################################### 
# USAGE: 
#     annotateSV.R -vcf <input.vcf> -ref <ref.genome> -out_vcf <output.vcf>
# Arguments:
#   input.vcf (-vcf) - Input VCF file with breakend notation (BND).
#   ref.genome (-ref) - reference genome build 'hg38', 'hg19', etc.
#   output.vcf (-out_vcf) - Optional output filename (default: <input_filename>_annot.vcf)
#   output.bedpe (-out_bedpe) - Optional output filename (default: <input_filename>_annot.bedpe)
#   vaf.thres (-vaf_thres) - VAF threshold for bedpe output (default = 0)


library(R.utils)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(rtracklayer)
library(stringr)


# Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), 'CTX', # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), 'INV',
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, 'INS', # TODO: improve classification of complex events
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == '-'), 'DEL',
                                     'DUP')))))
};


# Read command line arguments
ARGS <- commandArgs(
  trailingOnly = FALSE,
  asValues = TRUE,
  defaults = list(ref = 'hg38', vaf_thres = 0)
  );


file <- ARGS$vcf;
if (!file.exists(file)) {
  stop(sprintf('Input file %s does not exist!', file))
  };
  
if (is.null(ARGS$out_vcf) & is.null(ARGS$out_bedpe)) {
  out_vcf <- sub('.vcf$', '_annot.vcf', file)
} else {
  out_vcf <- ARGS$out_vcf
  };

#for testing
# vcf <- readVcf('../sample_data/gridss_SV/BLGSP-71-06-00285-01A-01D--BLGSP-71-06-00285-99A-01D--matched.gridss_somatic_filtered.vcf', 'hg38')
vcf <- readVcf(file, ARGS$ref)

info(VariantAnnotation::header(vcf)) <-
  unique(as(rbind(as.data.frame(info(VariantAnnotation::header(vcf))), 
                  data.frame(row.names=c('SIMPLE_TYPE', 'VAF'),
                             Number=c('1', '1'),
                             Type=c('String', 'Float'),
                             Description=c('Simple event type annotation based purely on breakend position and orientation.',
                                           'Variant Allele Fraction')
                             )), 'DataFrame'))

gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)

# Add values to info field
info(vcf)$SIMPLE_TYPE <- 'BND'
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype

info(vcf)$SVLEN <- NA
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

info(vcf)$VAF <- ifelse((info(vcf)$SVLEN > 1000 | is.na(info(vcf)$SVLEN)),
                        info(vcf)$VF / (info(vcf)$VF + info(vcf)$REF + info(vcf)$REFPAIR), 
                        info(vcf)$VF / (info(vcf)$VF + info(vcf)$REF))

#rename ID and MATEID by prefixing with SVTYPE
vcf_renamed <- vcf
names(rowRanges(vcf_renamed)) <- paste0(info(vcf_renamed)$SIMPLE_TYPE, '_', names(rowRanges(vcf_renamed)))
#rownames(info(vcf)) <- paste0(info(vcf)$SIMPLE_TYPE, '_', rownames(info(vcf)))
#info(vcf[gr$sourceId])$MATEID <- paste0(info(vcf[gr$sourceId])$SIMPLE_TYPE, '_', info(vcf[gr$sourceId])$MATEID)

# write VCF file
writeVcf(vcf_renamed, out_vcf)

gr$VAF <- info(vcf[gr$sourceId])$VAF

if (!is.null(ARGS$out_bedpe)) {
  gr_renamed <- subset(gr, VAF >= ARGS$vaf_thres)
  
  if (length(gr_renamed) > 0) {
    seqlevels(gr_renamed) <- gsub('chr', '', seqlevels(gr_renamed))
    names(gr_renamed) <- paste0(info(vcf[gr$sourceId])$SIMPLE_TYPE, '_', names(gr_renamed))
    gr_renamed$sourceId <- paste0(info(vcf[gr$sourceId])$SIMPLE_TYPE, '_', gr_renamed$sourceId)
    gr_renamed$partner <- paste0(info(vcf[gr$sourceId])$SIMPLE_TYPE, '_', gr_renamed$partner)
  }
  rtracklayer::export(breakpointgr2pairs(gr_renamed), con=ARGS$out_bedpe)
} else {
  write.table()
}


