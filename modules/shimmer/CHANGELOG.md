# Changelog

All notable changes to the `shimmer` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-09-10

This release was authored by Jasper.

Shimmer requires paired samples (tumour vs normal) and a reference file. The main output files include a somatic_diff.vcf (or a somatic_diff.vs file) and a somatic_indels.vs. (.vs is a VarSifter file format to show variants). The VarSifter file format includes more descriptions (including he reference and alternate allele, the q-value, the genotype [0/0 normal and 0/1 for tumour], and depth of coverage)

- Initially tried to implement a rule to run shimmer on a per-chromosome basis before merging the vcfs, but the --region command is broken(?) or unreliable in shimmer. It doesn't run to completion for some reason.
- Ultimately returned it back to a one-time full shimmer run on the entire bam file.

Output:
Vs. file:
Index, chr, leftFlank, rightFlank, ref_allele, var_allele, muttype, coverage_norm, coverage_tum, ratio_ReadsOfaltBase_to_totalReads_normal, ratio_ReadsOfaltBase_to_totalReads_tumour, q_value 
