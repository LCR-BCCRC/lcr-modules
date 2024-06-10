# Changelog

All notable changes to the `cnvkit` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [1.0] 2022-11-26

Added in a metrics rule to track the quality of the segmentation. It will output the MAD, number of segments, and IQR. Higher MAD would indicate more noise. This will be used to compare with metrics of other CNV calling tools.

# [1.0] - 2022-02-17

Bundling does not scale well with large samples/new samples. Split up all the commands into their individual components.

Notably, an option for "new_normals" must be filled in. This option is used to generate a new panel of normals. If done, everything will be re-ran. If not, the same panel of normals will be used for all new samples.

Using a genome fasta, cnvkit will get a set of accessible regions. It will by default hard-mask the genome. I have also taken steps to remove chrG, chrJ, chrM to keep everything consistent.

Cnvkit will get the capture space if filled in in the metadata. If not, it will use the default exomes. From there, a series of target and antitarget regions will be defined.

Coverage at these regions will be calculated for all samples.

A panel of normal will be built from the normals in the group.

From there, every tumour sample will be normalized and segmented based on the panel of normals.

BAF is determined by calling SNPs from dbSNP and using that as a further filtering step. Threshold is used as default instead of clonal, since ploidy and purity is not already known/determined for each sample.

CNVs are further filtered based on confidence intervals.

From there, the output is used to generate scatter plots and chromosome diagrams. Sex is inferred. Gene level CN is determined (at least for confident ones) - based on the intersect of segment-level CNV and bin-level CNV. It also tries to determine potential breakpoints. Finally, it outputs a seg file.

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
