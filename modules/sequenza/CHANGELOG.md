# Changelog

All notable changes to the `sequenza` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4] - 2026-02-11

Additional changes not requiring a version bump, existing files will be overwritten instead. Authored by Sierra Gillis.

- update cnv2igv step to have the `--preserve_log_ratio` flag and use 1.5 version of lcr-script; config has an option to remove this flag
- update lcr-script versions: liftover --> 2.0, fill_segments --> 1.2
- updated liftover to expect 0-based input and to write output as 0-based as well, and to fix dropping of large unmapped segments
- update fill_segments step to drop non-canonical chrs
- update normalize_projections step to format `NA`s correctly in output
- better formatting of rules and configs i.e. removed the need to specify the output paths in the config and cleaned up prepare_projection input function
- cleaned up formatting to be consistent between all CNV calling modules
- update rule all to prevent ambiguous wildcards
- better resource control

## [1.4] - 2021-01-17

This release was authored by Kostiantyn Dreval.

- This update includes several new features, such as: cram support, filtering of variants that falls in ENCODE blacklisted regions, version check of oncopipe module, and handling of input and output symlinks.

## [1.3] - 2020-07-17

This release was authored by Bruno Grande.

- The `_sequenza_merge_seqz` rule was concatenating all input files instead of just the SEQZ files. This bug was fixed.


## [1.2] - 2020-07-16

This release was authored by Bruno Grande.

- All `rules` references were wrapped with `str()`.

## [1.1] - 2020-07-16

This release was authored by Bruno Grande.

- I moved the GC content reference file from the `_sequenza_all` rule to earlier in the workflow because using `reference_files()` doesn't seem to work when you refer to a rule's input using `rules._sequenza_all.input` (as it's done in the demo project Snakefile).

## [1.0] - 2020-06-06

This release was authored by Ryan Morin and Bruno Grande.

<!-- TODO: Explain each important module design decision below. -->

- No module design decisions explained here yet.
