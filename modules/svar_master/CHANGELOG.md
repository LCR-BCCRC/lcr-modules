# Changelog

All notable changes to the `svar_master` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-12-29

This release was authored by Laura Hilton.

- This module integrates Manta, GRIDSS, and the HMFtools PURPLE-LINX tools into a single pipeline. It combines the Manta and GRIDSS SV calls into a single bedpe file and annotates the breakpoints by the nearest gene/locus based on a user-provided bed file. 
- Since this module includes the GRIDSS and Manta Snakefiles but also depends on the outputs of SLMS-3, SLMS-3 MUST be run before this module. It is also highly recommended to run GRIDSS on the majority of samples prior to launching the complete pipeline. This can be accomplished by specifying `_gridss_all` as the target when launching Snakemake.  
- Since this pipeline contains many, many steps it is also highly recommended to throttle resources carefully when running on a cluster. See the configs for each module for how resource-intensive jobs can be throttled. 


## [1.0] - 2023-11-06
This change was authored by Sierra Gillis.
Small update to the rule that creates bed files for annotation. It will only create bedfiles and annotate with bedtools when the input bedpe contains variants. Otherwise it creates empty bedfiles. The following rule that combines annotated bed files will also output an empty bedpe if the input bedfiles are empty. This will avoid the workflow throwing an error when samples to do have variants.