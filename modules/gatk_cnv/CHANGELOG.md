# Changelog

All notable changes to the `gatk_cnv` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [Unreleased]

## [1.0] (draft)

This release was authored by Ryan Morin.

- GATK4 somatic copy number workflow for exome (capture) data: PreprocessIntervals (padded targets, chromosomes matching `chromosome_regex`, autosomes by default), AnnotateIntervals, CollectReadCounts, one CreateReadCountPanelOfNormals per seq_type--genome_build--capture_space from all normals in the samples table except `pon_exclude`, CollectAllelicCounts at common gnomAD SNPs, DenoiseReadCounts, ModelSegments and CallCopyRatioSegments
- every GATK step is given the reference FASTA, so CRAM input works
- targets with denoised log2 copy ratio below `min_log2_copy_ratio` (-5) or above `max_log2_copy_ratio` (5, CN 64) are dropped before segmentation; the zero-coverage dropouts (common in FFPE) and high-side single-target artefacts otherwise each become a segment, the latter with CN in the hundreds
- ModelSegments defaults (changepoint penalty 5, one smoothing iteration per fit) come from tests on Chapuy FFPE exomes, where GATK's defaults gave 1,400-4,900 segments per tumour and these gave 31-121, with better agreement with the DLBClass ABSOLUTE segments
- paired tumours use heterozygous sites from the matched normal; unpaired tumours are run tumour-only (`no_normal`, sites genotyped from the tumour, minimum depth 30). Copy ratio is equally reliable either way; tumour-only allele fractions are not, in LOH regions of high-purity tumours
- optional genome-wide plot per tumour (`output: plot`), drawn by `src/plot_modelsegments.R` in an R env (the pinned PureCN env by default) rather than GATK's PlotModeledSegments, which needs R inside the GATK env; targets removed before segmentation are shown in grey
- standardized seg (ID chrom start end module LOH_flag CN log.ratio) with CN = round(2 * 2^log2) (not purity-corrected) and LOH_flag NA, lifted to grch37/hg38 projections, filled and normalized as in the other CNV modules
