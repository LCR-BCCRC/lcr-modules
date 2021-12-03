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

## [1.2] - 2020-09-10

This release was authored by Jasper.

Fixed a config bug (snakemake replacing 3&4 with 3numdegree4 in config_WGS.txt) - added \ on top of &. Note: tested it, the output was the same (since 3&4 is default)
Changed output dir to be controlfreec-1.1/ instead of controlfreec-1.0/
Added a new feature in snakemake to allow for conversion of CNV calls to bed.

## [1.2] - 2020-11-05

This release was authored by Jasper.

Improved control-freec calls for FFPE genomes. Everything is now default in paired mode with BAF turned on to measure LOH/BAF. All unmatched samples will be ran with a default normal. This will now generate temporary pileups in a scratch space (time consuming/rate limiting step). The parameters have been set to minimize noise while also maximizing the number of prominent signals across the genome. To accomplish this, I set the default ploidy to 2 to prevent control-freec from overfitting the model for unpaired samples when ploidy != 2. Importantly, I set the parameter forceGCcontentNormalization=1, which forces GC normalzation for both the sample and control read counts (RC), and then calculates the ratio of sample RC/control RC. This effectively smoothed out the profiles without added segments of overamplified "noise". I also allowed control-freec to assess the relative normal contamination and allowed control-freec to set its own window size given the coverage across the genome. These parameters effectively smoothed out unnecessary noise even further.

Notably, in paired mode, with BAF mode on, FREEC normalizes with GC-content, and without BAF mode on, FREEC normalizes using a control sample. Previously, I would get fairly clean profiles, with some marked "artifacts" spanning large chunks of the genome. These were ascribed as "artifacts" because I saw them in controlled samples as well. This may be a contribution of multiple things, including this spike-in used in library preparation we observed at the start of chr8 that mapped there and the uneven distribution of FFPE samples. I tried to implement noisy=T and increased minimum coverage and minimum CNA length to 8 to smooth out blips in signal, but those large artifacts remained. I tried to run these samples in unpaired mode to allow it to normalize itself with just GC content, but the artifacts still remained. Finally, I used the current implementation to decrease "noise" by normalizing samples by their own GC-content and then I let the samples normalize against each other to remove these excessive artifacts that appear in both samples vs controls.

This implementation has been tested on unmatched samples too using a high coverage, normal FFPE sample, and it has shown to display clean profiles in these cases too.

Note: this version is not meant for capture/exome data.

## [1.2] patch 2021-02-25
Added GEM mappability features - can now use/generate a hard-masked mappability file (useful for FFPE genomes) with the setting "hard_masked" = True. If this is set, GEM will be installed and ran on your reference genome of choice.

Also added freec2circos function.