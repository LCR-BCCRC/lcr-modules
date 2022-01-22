# Changelog

All notable changes to the `cnvkit` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-01-22

This release was authored by Jasper Wong.

Cnvkit uses both on-target and off-target bait regions to determine CNV.

Everything is bundled under one function "batch", which builds a normal reference from a panel of normals and calls CN for each tumour.

Copy number reference profile (.cnn) - reference files, tumour samples
- contains GC content, repeat-masker masked proportion of sequence region, statistical spread of dispersion
- for tumour samples, it will show the coverage depth of sample across target sites and anti-target sites

Bin-level log2 ratios (.cnr) - tumour samples
- chr, start, end, gene, log2 ratio, depth, weight 
- weight is generated from the deviation of the log2 ratio of the normal panel from the neutral coverage (i.e. from 0.0)

Segmented log2 ratios (.cns) - tumour samples
- same as .cnr, but bins are stitched together to form segments
- additional column probes, which indicates the number of bins covered by a segment

Segmented log2 ratios (.call.cns) - tumour samples
- same as .cns, but contains absolute copy number (based on log.ratio and a t-test)
- also contains a t-test p-value (notably it is a test against log.ratio 0.0) - if it is n.s., then it is assumed to be CN = 2

BAF was also added by using bcftools to call SNPs from a reference dbSNP file. This was then used to determine BAF.