# Changelog

All notable changes to the `hmftools` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-07-29

This release was authored by Laura Hilton.
- This module takes output from the `gridss` module to perform purity/ploidy/CNV calling (PURPLE) and generate SV plots (LINX). 
- Is able to run in unmatched normal mode, but the results will be substantially noisier than on matched samples. 
- For optimal pipeline efficiency, include both `modules/gridss/1.0/gridss.smk` and `modules/hmftools/1.0.smk` in the same Snakefile, with `_hmftools_all` as the target rule. 
