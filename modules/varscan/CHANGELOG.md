# Changelog

All notable changes to the `varscan` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-22

This release was authored by Helena Winata.


- The `varscan` module supports somatic variant calling for paired tumour-normal samples and germline variant calling for unpaired (`no_normal`) samples.
- Command line parameters for somatic and germline calls are specified for each `seq-type`.
- Intermediate mpileup files are temporary due to file size.
- Requires the `vcf2maf` module to convert VCF files to MAF files (final outputs).

