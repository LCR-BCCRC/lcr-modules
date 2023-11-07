# Changelog

All notable changes to the `spechla` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-10-26

This release was authored by Jasper Wong.

Beta version of running SpecsHLA - extracts reads from bam/cram from HLA locus and then deconvolutes them.
Compatible with capture, WGS, RNA-seq, long-read

- Currently works with WGS, WES, mRNA. Can run on long-reads too but just modify the settings.

The main caveat here is that the tool doesn't distinguish samples that are NULL because it loses the HLA or if it is NULL because it just cannot deconvolute the reads.