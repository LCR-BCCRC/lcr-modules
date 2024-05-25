# Changelog

All notable changes to the `combine_cnv` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-04-21

This release was authored by Kostiantyn Dreval. The first version of the cnv merging snakefile. It will generate a mster seg for a projection requested by user, and a parallel `contents` file containing paths to all individual seg files that went into the merge.

It can handle mixed sample tables with one or both of `genome`, `capture` samples present.

The individual CNV callers are specified by user in config under `names` key. The order of callers specified in `names` will also determine the preference of choosing the seg file for merge if outputs from more than one tool are available. If the user desires to merge outputs of only one tool, then only that specific caller should be specified under `names` key.

When this module is combined in project snakefile with individual CNV callers, it will trigger run of individual CNV caller module if the outputs hasn't been previously generated.
