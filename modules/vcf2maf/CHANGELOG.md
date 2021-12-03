# Changelog

All notable changes to the `vcf2maf` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2] - 2020-12-01

This release was authored by Kostia Dreval

- This version of modules hamdles vcf2maf feature of specifying non-canoniical ENST IDs to override canonical selection. It can be specified in config as a path to txt file containing list of IDs. If no IDs to be provided, the `switches` should be left blank. Separate lists can be provided for different genome builds. If file with custom transcripts is specified but does not have any transcripts listed, it should contain at least new line character as vcf2maf checks for this file to be more than 0 b. In addition, decompressed vcf files in this version are marked as `temp()` to be deleted after conversion, since they were left in the module folders and taking unnecessary disk space. Finally, resources restriction was enabled in module configuration, because multiple jobs using the same vep file created I/O bottleneck and slowed down some systems.


## [1.1] - 2020-11-04

This release was authored by Kostia Dreval

- 3rd level analyses of variants often utilize tools that are restricted to either hg19 or hg38 genome build. Therefore, this release of vcf2maf module includes automatic conversion between genome builds using crossmap tool. This implies that maf files in both original and converted genome builds are generated, and conveniently avaliable from the same output folder.


## [1.0] - 2020-06-05

This release was authored by Helena Winata.

- This module is designed as a utils-type module. It does not require `op.setup_module` as it is executed as part of a primary module that generates vcf files (varscan, manta, strelka, etc.)
- The config file is used to specify the vep_cache location and command line arguments.
