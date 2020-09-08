# Changelog

All notable changes to the `utils` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/2.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0] - 2020-05-02

This release was authored by Helena Winata.


- Usage of `CFG["dirs"]["parents"]` is removed as multiple includes of `utils.smk` are ignored and therefore, `CFG` are not updated.
- `CFG["dirs"]["parents"]` is replaced by `{out_dir}/{prefix}/{suffix}` where the wildcards are constrained such that `wildcards.prefix` always starts with 2 numbers (e.g. `02-bam_sort/path/to/file`)
- The variable `CONFIG` is used as a shortcut to `config["lcr-modules"]["utils"]` as the utils module does not rely on `op.setup_module`
- `wildcards.prefix` is also constrained withing each rule to ensure that the symlinking rules within the primary module (STAR or bwa_mem) are executed. It ensures that `utils.smk` rules only try to access inputs from a certain directory (e.g. sorting rule will only look in `02-bam_sort`)
- The snakefile does not include a target rule as it relies on the primary module to call final outputs.

