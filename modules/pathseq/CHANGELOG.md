# Changelog

All notable changes to the `pathseq` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-03-18

This release was authored by Kostiantyn Dreval.

- The reference files are downloaded directly from GATK where available. The grch37 bundle was discontinued and therefore I followed their recommendations to develop rules generating reference files in a consistent way. It also required to remove EBV chromosome from hg38 reference.

- The tool is very resource-consuming, therefore RAM requirements are very high.

- All conda environment files are reused from other modules to reduce the number of files stored and avoid duplication of environments containing same tools.

- The cut-off for assigning EBV status can be provided through config and will be used to determine EBV-positivity. It should be provided as a fraction of EBV reads relative to the total number of mapped reads.

- Currently it is designed to run on all samples, including bam/cram, rnaseq/genome, and tumour/normal. This can be modified in later versions as needed.
