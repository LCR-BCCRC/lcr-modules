# Changelog

All notable changes to the `whatshap` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - ???

This release was authored by Nicole Thomas.

## [2.0] - 2024-10-08

This release was authored by Laura Hilton. 

- This is a more generalized version of the WhatsHap module. 
- Introduces parallelization of the phase_vcf rule across chromosomes. 
- Enables the use of a VCF derived from matched short-read sequencing alongside long-read alignments, matchin based on biopsy_id. 
- Can take a bed file of regions-of-interest to create a phased bam so that whole genome phased bams are optional. 
