# Changelog

All notable changes to the `dnds` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-08-02

This release was authored by Kostia Dreval.

- The version 1 uses the reference files supplied directly with dNdS.
- Currently, only the grch37-based maf files are supported.
- The snakefiles relies on the lcr-script `generate_smg_inputs` to produce the input maf.
- The implementation is similar to the other SMG lcr-modules, e.g. mutsig.
