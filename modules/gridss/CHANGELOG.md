# Changelog

All notable changes to the `gridss` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-07-22

This release was authored by Laura Hilton.
- Separates GRIDSS preprocessing steps per sample rather than per tumour-normal pair to speed up this step and eliminate redundant normal preprocessing. 
- Groups GRIDSS preprocessing with GRIDSS run steps to promote deletion of temp preprocessing files as quickly as possible. For this reason it is highly rceommended that this pipeline be run on a cluster. 
- IMPORTANT: If more than one tumour will be run with an unmatched normal, it is highly recommended to run the _gridss_preprocess step with the unmatched normal wildcards as the target BEFORE launching the rest of the pipeline on the tumours and tumour/normal pairs. A more sustainable fix to this will come in a future release. 
- Retains a no-normal pipeline for e.g. specific capture data. 
