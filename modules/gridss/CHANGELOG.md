# Changelog

All notable changes to the `gridss` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-07-22

This release was authored by Laura Hilton.
- Separates GRIDSS preprocessing steps by tumour and normal to speed up this step and eliminate redundant normal preprocessing. 
- Groups GRIDSS preprocessing with GRIDSS run steps to promote deletion of temp preprocessing files as quickly as possible. For this reason it is highly rceommended that this pipeline be run on a cluster. 
- Retains a no-normal pipeline for e.g. specific capture data. 
