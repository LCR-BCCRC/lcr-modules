# Changelog

All notable changes to the `qc` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-05-18

This release was authored by Kostia Dreval.

<!-- TODO: Explain each important module design decision below. -->

- This module replaces the current `picard_qc` module.
- Verion 1 uses same bed file for both baits and target regions. This can be separated in the future versions if this functiionality will be of interest.
- The module generates supplementary table-ready tab-deliminated file with QC metrics to be reported as a good practise described by MIRAGE.
- The `MeanCorrectedCoverage` means that the filters for unmapped, low-quality, duplicated, secondary, and supplemental reads has been applied before calculating coverage. The reported value is effective coverage.
- Additional flags for each tool can be specified via config.
