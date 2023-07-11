# Changelog

All notable changes to the `fishhook` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-06-27

This release was authored by Jacky Yiu.

- Both grch37-based and grch38-based maf files are supported, please supply coveriate file of the same build
- The snakefiles relies on the lcr-script `generate_smg_inputs` to produce the input maf.
- The implementation is similar to the other SMG lcr-modules, e.g. mutsig, dnds.
