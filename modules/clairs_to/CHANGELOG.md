# Changelog

All notable changes to the `clairs_to` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-02-20

This release was authored by Giuliano Banco.

This module is designed to only work with a sample table that includes the columns 'chemistry' and 'platform'. These are required by schemas. Note that ClairS-TO was developed to support R10 data but not R9, however, you can specify R10 as the chemistry value and use the R10 models.

This module only works with PromethION and hg38 data, as of version 1.0.