# Changelog

All notable changes to the `controlfreec` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-07-29

This release was authored by Jasper.

Basic controlFREEC:
-input: 
    - WGS or WES bam
-output: 
    - text files indicating CNV positions, general characteristics of file (ex. ploidy, sample purity)
    - graphs visualizing CNV positions per chromosome

- Implemented for use with WGS libraries.
- Compiled one environment containing samtools, sambamba, bedtools in freec.
- Features the most basic default control-freec run (generates CNV plots, CNV tables, and sublone txt).

- Config file:
- There are 2-5 parameters: [general], [sample], 
    - [general] - parameters used for controlfreec - there are different defaults for WGS vs. WES
    - [sample] - parameters for bam files used (also path/to/bam)
* version 1.0 only allows for unpaired mode so far

- (optional: [control], [BAF] and [target].)
    - [control] - path/to/control
    - [BAF] - for calculating BAF profiles and call genotypes (to detect LOH)
    - [target] - provide a .bed file with coordinates of probes, exons, amplicons for exome-sequencing or targeted-sequencing. Set "window=0" in [general] to use read count "per exon" instead of "per window"
Additional features in config can be found here (http://boevalab.inf.ethz.ch/FREEC/tutorial.html#CONFIG)

## [1.1] - 2020-09-04

This release was authored by Jasper.

Updated controlfreec.smk to allow for running on Numbers cluster.

Fixed bugs related to genome builds and preparation of controlfreec-specific reference files. Now will generate chr.len file for hg19, hg38 (w/ chr prefix), grch37, grch38 (no chr).