# Changelog

All notable changes to the `mosdepth` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-04-22

This release was authored by Giuliano Banco.

- This module uses Mosdepth to calculate sequencing coverage from input BAM/CRAM files

- This module supports optional calculation of per-base or region-based coverage depending on configuration, allowing users to determine coverage within targeted sequencing panels.

- This module supports optional suppression of determining per-base coverage using the `no_per_base` setting in the config. This reduces runtime when users do not which to identify per-base coverage.
