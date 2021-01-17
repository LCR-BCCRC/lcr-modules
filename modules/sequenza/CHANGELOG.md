# Changelog

All notable changes to the `sequenza` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.5] - 2021-01-17

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
