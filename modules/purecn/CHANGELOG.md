# Changelog

All notable changes to the `purecn` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-07-06

This release was authored by Jasper.

- Built the module with the intention of integrating it with cnvkit.
- PSCBS is the default segmentation method, and leverages variant information in conjunction with the CBS algorithm to call more accurate CNVs.
- PureCN uses both germline variants (from normals) and somatic variants to have more confidence in calling CNVs and to also filter out noise.

The overall workflow of the module is as follows:
Step 1: set up references
- This module is built with the intention of using after cnvkit was ran on capture libraries. The module first sets up reference files, which includes symlinking CNVkit coverage files and GEM mappability files.
- PureCN then sets up interval files based on the specific capture space of the capture kit. I customized the intervalFile.R script (in etc/) to not require the fasta index to be in the base path of the fasta file (i.e. we symlink it in our reference workflow).

Step 2: set up normals
- PureCN then calls germline variants using Mutect2.
- A vcf containing a panel of normals is constructed (this would be used for both the database below and also for the normal database used by PureCN)
- From there, a panel of normals database is constructed using Mutect2 GenomicsDBImport and CreateSomaticPanelOfNormals (to be used by Mutect2)

Step 3: set up coverage files
- Coverage files were generated. This part was trickier since pureCN does not allow for CRAM files. Rsamtools have CRAM compatibility in newest versions but it requires R at least 4.2 and gphost08 was not.
- I tried to integrate bamUtils.R's version of CRAM importing by calling samtools in the command line within conda environment, but it got a bit too complicated in the code.
- In the end, I opted to use GATK's depthOfCoverage, which was already compatible with pureCN
- PureCN coverage.R was customized in etc/ to allow for sample_id calling as input and output.
- PureCN's coverage.R script was ultimately just used to GC normalize the GATK depthOfCoverage coverage files.

Step 4: call variants in tumors
- Mutect2 was ran on all tumor samples, with the panel of normals (from the database above) used as a way to annotate whether or not this germline variant was observed in the tumor samples. Ultimately, the variant file that comes out will have both somatic and germline mutations annotated (which is why we couldn't do a simple mutect2 unpaired mode run). PureCN requires this specific format.

Step 5: call PureCN
- This part is bifurcated into two streams:

1. PureCN using coverage from CNVkit but just with the PSCBS algorithm. (PureCN_cnvkit)

2. PureCN using coverage called and normalized by its interal scripts (PureCN_denovo)

- in both cases, a mapping bias database is constructed using a panel of normals

1. PureCN_cnvkit just uses the panel of normal vcfs to create mapping bias database

2. PureCN_denovo uses both panel of normal vcfs and the coverages (from GATK/pureCN) to create a mapping bias and a normalDB database

- an additional tidy shell script (in etc/) was used to set up the seg files for downstream processes

Step 6: Project to other genome spaces
- The last part is to project the seg files into other genome builds and to also remove the wildcard "capture_space" to keep things consistent with LCR-module wildcard conventions.

- The reason there's a split is because some samples seem to work better with CNVkit + PureCN PSCBS and some samples work better alone with just pureCN PSCBS (allowing pureCN to segment it).
- An additional step (beyond this module) would be to select the best seg file (i.e. the least noisy seg file) out of the three exome CNV callers: (cnvkit, pureCN with cnvkit coverage and segmentation, and pureCN only)
- This would require comparing MAD (median absolute deviations) scores - the lowest score wins.
