#!/usr/bin/env Rscript
# Wrapper around PureCN::runAbsoluteCN() for GATK4 CNV outputs.
#
# Configuration matches DepMap's PureCN.wdl invocation (broadinstitute/depmap_omics):
#   * Mutect2 VCF with --genotype-germline-sites=true (BAF anchors for
#     segmentation and purity refinement)
#   * IntervalFile.R-generated intervals (Gene + gc_bias + mappability columns)
#     so callAlterations() produces real per-gene absolute CN
#   * Hclust segmentation, post.optimize=TRUE, model.homozygous=TRUE
#   * Cell-line-tight purity bounds (default 0.90 - 0.99)
#
# Inputs (via --flags):
#   --seg_modelFinal  GATK ModelSegments output (.modelFinal.seg)
#   --mutect_vcf      Mutect2 VCF (with FilterMutectCalls applied);
#                     empty/zero-byte file → coverage-only mode
#   --intervals       PureCN-format interval file (from IntervalFile.R)
#   --sample_id       Sample identifier (used as output prefix)
#   --out_dir         Output directory
#   --genome          "hg19" / "hg38" (PureCN's expected labels)
#   --fun_segmentation  Hclust / PSCBS / CBS / GATK4
#   --min_purity, --max_purity, --max_copy_number, --max_segments  bounds

suppressPackageStartupMessages({
  library(optparse)
  library(PureCN)
})

option_list <- list(
  make_option("--seg_modelFinal",  type = "character"),
  make_option("--mutect_vcf",      type = "character", default = ""),
  make_option("--intervals",       type = "character"),
  make_option("--sample_id",       type = "character"),
  make_option("--out_dir",         type = "character"),
  make_option("--genome",          type = "character", default = "hg19"),
  make_option("--fun_segmentation", type = "character", default = "Hclust"),
  make_option("--min_purity",      type = "numeric", default = 0.90),
  make_option("--max_purity",      type = "numeric", default = 0.99),
  make_option("--max_copy_number", type = "integer", default = 8),
  make_option("--max_segments",    type = "integer", default = 1000)
)

opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Treat zero-byte placeholders as "no VCF" so coverage-only mode still works
# if the Mutect2 step is skipped or fails for a sample.
has_vcf <- nzchar(opt$mutect_vcf) &&
  file.exists(opt$mutect_vcf) &&
  file.info(opt$mutect_vcf)$size > 0

if (!has_vcf && opt$fun_segmentation == "PSCBS") {
  message("No Mutect2 VCF supplied; falling back to fun.segmentation='Hclust'")
  opt$fun_segmentation <- "Hclust"
}

# runAbsoluteCN's fun.segmentation parameter expects a function reference.
fun_segmentation <- switch(
  opt$fun_segmentation,
  "CBS"    = PureCN::segmentationCBS,
  "PSCBS"  = PureCN::segmentationPSCBS,
  "Hclust" = PureCN::segmentationHclust,
  "GATK4"  = PureCN::segmentationGATK4,
  stop("Unknown fun_segmentation: ", opt$fun_segmentation)
)

args <- list(
  seg.file         = opt$seg_modelFinal,
  interval.file    = opt$intervals,
  sampleid         = opt$sample_id,
  genome           = opt$genome,
  fun.segmentation = fun_segmentation,
  test.purity      = seq(opt$min_purity, opt$max_purity, by = 0.01),
  test.num.copy    = 0:opt$max_copy_number,
  max.segments     = opt$max_segments,
  # DepMap's PureCN config — matches `--post-optimize --model-homozygous`.
  # post.optimize re-fits purity/ploidy after initial grid search.
  # model.homozygous treats LOH segments as truly homozygous (no minor allele),
  # which fits clonal cell-line genetics better than mixed-population assumptions.
  post.optimize    = TRUE,
  model.homozygous = TRUE,
  plot.cnv         = FALSE
)

if (has_vcf) {
  args$vcf.file <- opt$mutect_vcf
}

ret <- tryCatch(
  do.call(runAbsoluteCN, args),
  error = function(e) {
    message("runAbsoluteCN failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

file_prefix <- file.path(opt$out_dir, opt$sample_id)

# Purity / ploidy summary
purity_df <- data.frame(
  sample_id = opt$sample_id,
  purity    = ret$results[[1]]$purity,
  ploidy    = ret$results[[1]]$ploidy,
  flagged   = if (!is.null(ret$results[[1]]$flagged)) ret$results[[1]]$flagged else NA,
  curated   = if (!is.null(ret$results[[1]]$curated)) ret$results[[1]]$curated else NA,
  stringsAsFactors = FALSE
)
write.csv(purity_df, paste0(file_prefix, ".csv"), row.names = FALSE)

# Per-gene absolute CN — works now that the interval file has a Gene column.
genes_df <- tryCatch(
  callAlterations(ret),
  error = function(e) NULL
)
if (!is.null(genes_df)) {
  write.csv(genes_df, paste0(file_prefix, "_genes.csv"), row.names = TRUE)
} else {
  write.csv(
    data.frame(gene = character(), C = integer(), seg.mean = numeric(), seg.id = integer()),
    paste0(file_prefix, "_genes.csv"),
    row.names = FALSE
  )
}

loh_df <- tryCatch(
  callLOH(ret),
  error = function(e) NULL
)
if (!is.null(loh_df)) {
  write.csv(as.data.frame(loh_df), paste0(file_prefix, "_loh.csv"), row.names = FALSE)
} else {
  write.csv(
    data.frame(chr = character(), start = integer(), end = integer(), type = character()),
    paste0(file_prefix, "_loh.csv"),
    row.names = FALSE
  )
}

seg_df <- tryCatch(
  predictSomatic(ret),
  error = function(e) NULL
)
if (!is.null(seg_df)) {
  write.csv(as.data.frame(seg_df), paste0(file_prefix, "_segmentation.csv"), row.names = FALSE)
}

message("PureCN run complete: ", opt$sample_id,
        " purity=", round(purity_df$purity, 3),
        " ploidy=", round(purity_df$ploidy, 3))
