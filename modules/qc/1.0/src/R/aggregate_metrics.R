#!/usr/bin/env Rscript
#

suppressWarnings(suppressPackageStartupMessages({

  message("Loading packages...")

  required_packages <- c(
    "data.table",
    "tidyverse"
  )

  for (pkg in required_packages){
    library(pkg, character.only = TRUE)
  }
}))



# Collect the base quality
base_scores = fread(
  snakemake@input[[2]],
  skip = 7
) %>%
  as.data.frame

QualityScore =
  base_scores %>%
  mutate(
    totalQ = as.numeric(QUALITY) * as.numeric(COUNT_OF_Q),
    totalQ = sum(totalQ),
    totalN = sum(COUNT_OF_Q)
  ) %>%
  mutate(AverageBaseQuality = round(totalQ / totalN, 2)) %>%
  distinct(AverageBaseQuality) %>%
  `rownames<-`(snakemake@wildcards$sample_id)


# handle samtools metrics
command = paste("grep ^SN", snakemake@input[[1]], "| cut -f 2-")
samtools_scores = fread(cmd = command,fill=TRUE,sep = "\t") %>%
  select(-V3) %>%
  mutate(V1=gsub(":", "", V1)) %>%
  column_to_rownames("V1") %>%
  `names<-`(snakemake@wildcards$sample_id) %>%
  t %>%
  as.data.frame

SamtoolsMetrics =
  samtools_scores %>%
  select(
    `insert size average`,
    `average length`,
    `pairs on different chromosomes`,
    `raw total sequences`,
    `reads mapped`,
    `reads unmapped`,
    `reads duplicated`
  )

SamtoolsMetrics =
  SamtoolsMetrics %>%
  mutate(
    ProportionReadsDuplicated = round(`reads duplicated` / `raw total sequences`, 2),
    ProportionReadsMapped = round(`reads mapped` / `raw total sequences`, 2)
  )

colnames(SamtoolsMetrics) = c(
  "AverageInsertSize",
  "AverageReadLength",
  "PairsOnDiffCHR",
  "TotalReads",
  "TotalUniquelyMapped",
  "TotalUnmappedreads",
  "TotalDuplicatedreads",
  "ProportionReadsDuplicated",
  "ProportionReadsMapped"
)

# handle coverage metrics
if (snakemake@wildcards$seq_type == "capture"){
  command = paste("grep BAIT_SET", snakemake@input[[3]] ,"-A 1")
} else {
  command = paste("grep GENOME_TERRITORY", snakemake@input[[3]] , "-A 1")
}

CoverageMetrics = fread(cmd=command) %>%
  as.data.frame %>%
  `rownames<-`(snakemake@wildcards$sample_id)

if (snakemake@wildcards$seq_type == "capture"){
  CoverageMetrics =
    CoverageMetrics %>%
    select(
      BAIT_SET,
      BAIT_TERRITORY,
      MEAN_TARGET_COVERAGE,
      ZERO_CVG_TARGETS_PCT,
      PCT_TARGET_BASES_10X,
      PCT_TARGET_BASES_30X,
      FOLD_ENRICHMENT
    )
} else {
  CoverageMetrics =
    CoverageMetrics %>%
    select(
      GENOME_TERRITORY,
      MEAN_COVERAGE,
      PCT_10X,
      PCT_30X
    ) %>%
    mutate(
      ProportionTargetsNoCoverage=NA, .after = MEAN_COVERAGE,
      BAIT_SET = "whole-genome", .before = GENOME_TERRITORY,
      FOLD_ENRICHMENT = NA, .after = PCT_30X
    )
}

colnames(CoverageMetrics) = c(
  "TargetSpace",
  "Target_Territory",
  "MeanCorrectedCoverage",
  "ProportionTargetsNoCoverage",
  "ProportionCoverage10x",
  "ProportionCoverage30x",
  "FoldEnrichment"
)

print(colnames(CoverageMetrics))

outputMetrics = cbind(SamtoolsMetrics, CoverageMetrics, QualityScore) %>%
  relocate(AverageBaseQuality)

outputMetrics$SeqType = snakemake@wildcards$seq_type

outputMetrics =
  outputMetrics %>%
  relocate(SeqType) %>%
  as.data.frame %>%
  rownames_to_column("UID")

write_tsv(outputMetrics, snakemake@output[[1]])
