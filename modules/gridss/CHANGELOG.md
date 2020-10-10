# Changelog

All notable changes to the `gridss` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2020-10-09
This release was authored by Laura Hilton. See the [GRIDSS man page](https://github.com/PapenfussLab/gridss) for extensive documentation. 
- Add automatic reference file downloading from files hosted at the BCGSC [downloads page](https://bcgsc.ca/downloads/morinlab/hmftools-references/gridss/).
- Updated to hmftools gripss v 1.8.  


## [1.0] - 2020-07-22

This release was authored by Laura Hilton. See the [GRIDSS man page](https://github.com/PapenfussLab/gridss) for extensive documentation. 
- Separates GRIDSS preprocessing steps per sample rather than per tumour-normal pair to speed up this step and eliminate redundant normal preprocessing. 
- Groups GRIDSS preprocessing with GRIDSS run steps to promote deletion of temp preprocessing files as quickly as possible. However, job grouping only works on clusters so it is highly recommended that this pipeline be run on a cluster. 
- Retains a no-normal pipeline for e.g. specific capture data. The output of the no-normal pipeline is the complete SV vcf file, annotated with any viral sequence identified, and filtered for variants that pass GRIDSS filters: `gridss_viral_annotation_filtered.vcf.gz`.
