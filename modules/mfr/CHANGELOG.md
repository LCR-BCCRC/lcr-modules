# Changelog

All notable changes to the `mfr` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

[1.0] - 2026-09-14

This release was authored by Giuliano Banco.

Initial implementation of the mfr module (renamed from the mutation_foci prototype, with the version reset to 1.0 as a fresh identity). The module clusters non-coding mutation positions into "foci" via hierarchical clustering, scattered by sample_set and chromosome, selecting the cut height that maximizes mean silhouette width.

Inputs are per-sample, bgzip + tabix-indexed MAFs (inputs.sample_maf) rather than a single genome-wide master MAF, so a sample_set's jobs only ever open the MAFs of samples actually in that set. A new _mfr_extract_chrom rule (src/python/extract_chrom_maf.py) streams each sample's MAF through tabix {sample}.maf.gz {chrom} one line at a time, dropping coding Variant_Classification rows inline, so a chromosome's rows are never materialized as a table. cluster_foci.R then operates on that already-scoped file, no longer filtering by chromosome itself.

cluster_foci.R splits each chromosome's unique positions into gap-delimited chunks (wherever two consecutive positions are more than h_max apart) and clusters each chunk independently, rather than running one dist()/hclust() over every position. This is exact for any cut height <= h_max and bounds the O(n^2) dist() cost to the largest chunk instead of the whole chromosome, which matters when a sample_set pools ~700 WGS samples onto one chromosome.

A single combined conda env / container image (mfr: R clustering deps + Python + htslib) is used by both the extraction and clustering rules.