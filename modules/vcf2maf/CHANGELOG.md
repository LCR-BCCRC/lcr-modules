# Changelog

All notable changes to the `vcf2maf` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-05

This release was authored by Helena Winata.

- This module is designed as a utils-type module. It does not require `op.setup_module` as it is executed as part of a primary module that generates vcf files (varscan, manta, strelka, etc.)
- The config file is used to specify the vep_cache location and command line arguments.
