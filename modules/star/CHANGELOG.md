# Changelog

All notable changes to the `star` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-04-25

This release was authored by Bruno Grande.

<!-- TODO: Explain each important module design decision below. -->

- `utils.smk`
- `os.remove` is used instead of `temp()` since we can't modify the `utils` rules and we want to make sure the new BAM exists before deleting the old BAM
  the `utils` rules
- `new_bam` is there to ensure that the `old_bam` is deleted only once the `new_bam` exists
- chimeric reads
