# Changelog

All notable changes to the `lymphgen` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-11-02

This release was authored by Chris "coolbeans" Rushton.

Initial release, adding LymphGen, LGenIC, and allowing for CNV data (if availible)
You can run LymphGen using just SNVs, or with CNV and SV data
Note if you provide CNV and SV files, you should specify the appropriate column names in the config file
All possible iterations of LymphGen will be run (i.e. if you provide both CNVs and SNVs, LymphGen will be run
with both CNVs and SNVs, as well as just with SNVs)

## [1.0] - 2022-05-10

Additional improvements, authored by Chris "scienceparrot" Rushton.

Add an additional iteration of LymphGen, which excludes the A53 subgroup from classification


## [1.1] - 2024-08-23

This release was authored by Laura Hilton. 

- Significant update that uses as user-provided gene list for LymphGen instead of inferring target space from the input maf file. 
- Installs an updated version of [LGenIC](https://github.com/LCR-BCCRC/LGenIC/releases/tag/2.0.1) with a tagged release 2.0.1. 

## [2.0] - 2023-06-19

This release was authored by Faraneh Moayyed.

Runs LymphGen on a per-sample basis, allowing the user to provide individual maf files per sample. 
*Not recommended to be run in targeted mode because target space is inferred from the input maf file. This will create unexpected variation in the model across samples. Will also error out on empty maf files.*

## [2.1] - 2024-08-23

This release was authored by Laura Hilton. 

Significant updates include changes to [LGenIC](https://github.com/LCR-BCCRC/LGenIC/releases/tag/2.0.1) to accept an input gene list for targeted sequencing experiments, providing consistent results across samples. Minor modifications to the LymphGen run script ensure that the same target gene list is used for both CNV and SNV models, and allows samples with empty maf files to pass without error. 