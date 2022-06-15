# Changelog

All notable changes to the `mutsig` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-06-13

This release was authored by Kostia Dreval.

- MutSig2CV module relies on the "master maf" file containing the mutation data for all samples,
and a "sample_set" file where first column (sample_id or Tumor_Sample_Barcode) contains unique sample IDs,
and each column sets the name of the sample subset to be included in the run. The column for each
sample subset contains 1 if sample is a part of subset, or 0 if the sample is not part of the subset.

- Multiple input master maf files are allowed, for example when WGS and WES data are combined together.

- MutSig2CV is picky about the running directory and where all reference files are expected to be without
allowing user to modify or specify their location. This module will execute all runs in the `01-mcr` directory
and is configured to follow the standard lcr-modules logich of output files, so the results for each sample subset
are located in the `02-mutsig` directory and the main output is symlinked to `99-outputs`.

- MutSig2CV saves maf as part of it's output. Due to space concerns and because maf files can easily become massive
for large data sets, the maf files are meant to be temporary and are deleted upon completion of the module.

- If the particular sample subset is modified to add or exclude samples, please create a new column in the
sample_subsets file rather than modifying existing.
