# Changelog

All notable changes to the `purecn` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# [2.0] 2026-06-04

This release was authored by Sierra Gillis.

- panel of normals files are listed in the config and symlinked in this module from the panel_of_normals module
- runs `modes` using a wildcard instead of separate rules for each
- the `run` step is now a checkpoint that will try the command with different segmentation functions, and the the one that runs successfully will be symlinked for the downstream
- updated cnv2igv step to have the `--preserve_log_ratio` flag and use 1.5 version of lcr-script; config has an option to remove this flag
- updated lcr-script versions: liftover --> 2.0, fill_segments --> 1.2
- updated liftover to expect 0-based input and to write output as 0-based as well, and to fix dropping of large unmapped segments
- updated fill_segments step to drop non-canonical chrs
- updated normalize_projections step to format `NA`s correctly in output
- better formatting of rules and configs i.e. removed the need to specify the output paths in the config and cleaned up prepare_projection input function
- cleaned up formatting to be consistent between all CNV calling modules
- updated rule all to prevent ambiguous wildcards

# [1.0] Never fully PR'd

This version is kept here for legacy reasons, and some results have been generated with this version despite it not being PR'd