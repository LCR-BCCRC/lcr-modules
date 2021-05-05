# Changelog

All notable changes to the `ichorcna` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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