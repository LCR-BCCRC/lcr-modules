# Changelog

All notable changes to the `bwa_mem` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-04

This release was authored by Helena Winata.


- The `bwa_mem` module is designed to work with the `utils` module.
- `bwa_mem.smk` contain multiple symlinking rules to move `utils.smk` output to the next directory (e.g. sorted bam files are moved to `mark_dups` to be used as input to a rule that mark duplicates in `utils.smk`)
- Final output is a bam file with a naming format: `{sample_id}.{suffix}.bam` where the suffix is specified in `CFG["options"]["suffix"]`. The final output is symlinked to `99-outputs`.


## [1.1] - 2020-10-30

This release was authored by Laura Hilton.


- Included the `utils/2.0/utils.smk` module by default. 
- Resource unpacking is included for flexibility in specifying resource requirements. 