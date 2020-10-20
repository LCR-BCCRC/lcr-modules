# Changelog

All notable changes to the `star` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4] - 2020-10-20

This release was authored by Laura Hilton. 

- Read length is obtained from the samples table so that STAR can vary the index and `sjdboverhang` parameters per sample. The samples table MUST contain a `read_length` column. If all samples have the same read length, a default value can be set by following the instructions in the [LCR modules documentation](https://lcr-modules.readthedocs.io/en/latest/for_users.html#adding-and-transforming-columns).

## [1.3] - 2020-08-20

This release was authored by Kostiantyn Dreval.

- Rule `_star_run` was exiting with error on some systems when trying to remove temporary folder. `rmdir` from this rule was removed to prevent this error. 

## [1.2] - 2020-07-16

This release was authored by Bruno Grande.

- All `rules` references were wrapped with `str()`.

## [1.1] - 2020-06-06

This release was authored by Bruno Grande.

- The default configuration file was updated to include '__UPDATE__'.

## [1.0] - 2020-04-25

This release was authored by Bruno Grande.

- The BAM processing rules that would be generally useful were stored in `utils.smk`, which is a module intended to be shared between modules. 
- `os.remove` is used instead of `temp()` since we can't modify the `utils` rules and we want to make sure the new BAM exists before deleting the old BAM
  the `utils` rules. 
- When the symlinking happens, the previous version of the BAM file gets deleted to minimize the space usage. 
- By default, the chimeric reads are included in the BAM output files as soft-clipped reads (unlike the STAR default, which stores them in a separate file). See `--chimOutType WithinBAM SoftClip`. 
- By default, the original output STAR BAM files and sorted BAM files will be stored in scratch space to avoid clogging snapshots. 
