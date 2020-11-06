# Changelog

All notable changes to the `vcf2maf` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2] - 2020-11-04

This release was authored by Kostia Dreval

- 3rd level analyses of variants often utilize tools that are restricted to either hg19 or hg38 genome build. Therefore, this release of vcf2maf module includes automatic conversion between genome builds using crossmap tool. This implies that maf files in both original and converted genome builds are generated, and conveniently avaliable from the same output folder.


## [1.0] - 2020-06-05

This release was authored by Helena Winata.

- This module is designed as a utils-type module. It does not require `op.setup_module` as it is executed as part of a primary module that generates vcf files (varscan, manta, strelka, etc.)
- The config file is used to specify the vep_cache location and command line arguments.
