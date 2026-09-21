# Changelog

All notable changes to the `ichorcna` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2022-01-08

This release was authored by jawong.

IchorCNA has been updated to be compatible with .cram files.

Previously the issue with IchorCNA is the use of hmmcopy_utils ReadCounter, which is an old C tool that only works with bam files. With this version, I leverage deeptools BamCoverage to generate bigwig tracks from bam/cram, then I use UCSC bigWigToWig to lift it back to a wig file.

Notably, deeptool's bamCoverage stitches regions together that have identical values (i.e. the centromeres, which are all denoted as 0). This causes issues with ichorCNA since the reference .wig files have a specific number of lines (corresponding to the number of bins of X size across each reference genome). Also, the wig format used in ichorCNA requires a very specific notation (it is hard coded in their software). "fixedStep chrom=1 start=1 step=1000000 span=1000000" as a header. Every line after that is just the coverage value of that bin.

The UCSC BigWigToWig command converts the bigwig to essentially a bedGraph format (chr,start,end,coverage). Therefore, I leveraged bedOps --chop to slice the bedGraph into windows of binSize (wildcard), and then I used bedtools to intersect the bedGraph back to the sliced regions (thus maintaining the coverage value of the wig file, which gets deleted by bedOps).

Note: setting up ucsc-bigWigToWig conda env may not work for outdated OS (ex. numbers) 
(Problem: nothing provides __glibc >=2.17 needed by libgcc-ng-9.3.0-h5101ec6_17
Problem: nothing provides __glibc >=2.17 needed by libstdcxx-ng-9.3.0-hd4cf53a_17)
Run this on a server that is compatible first to set up the env (ex. gphost) then you can launch it on numbers.


## [1.0] - 2021-03-31

This release was authored by jawong.

IchorCNA requires indices to be ".bam.bai" format. Modified the initial symlink rule appropriately.
Modified runIchorCNA.R - line 128: seqinfo <- NULL - to allow for flexibility in genome builds (i.e. hs37d5).

IchorCNA is currently set up to run in unpaired mode using a panel of normals. It currently does not work with CRAM files.

The output files in 99-outputs include:
- .cna.seg - per-bin CNA state and log ratio
- .seg - segments called by the Viterbi algorithm - IGV compatible
- .segTxt - same as .seg but also includes subclonal status of segments (0 = clonal; 1 = subclonal) - not IGV compatible
- .genomeWide.pdf - CNA profile
- .param.txt - parameters including ploidy, tumour fraction, sex
- .corrDepth.txt - per-bin size corrected depth