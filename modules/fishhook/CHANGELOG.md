# Changelog

All notable changes to the `fishhook` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2023-06-27

This release was authored by Jacky Yiu.

- Both grch37-based and grch38-based maf files are supported, please supply covariates file of the same build
- The snakefiles relies on the lcr-script `generate_smg_inputs` to produce the input maf.
- The implementation is similar to the other SMG lcr-modules, e.g. mutsig, dnds.

## [1.1] - 2024-04-12
This update was authored by Sierra Gillis.

- Preprocessing maf step was updated to be agnostic of the input metadata samples dataframe. The values to subset this metadata by in order to get the IDs of samples in the sample set will be done through the level3_subsetting_categories.tsv table.
- The sample set outputs are tracked by launch-date of the snakemake and md5sum of the sample ids as well as sample set name. This implementation follows that in gistic2/1.1/

## [1.1] - 2024-06-26
This update was authored by Sierra Gillis.

- The "gene list" option for the model now takes the `gencode-33/gencode.annotation.grch37.gtf` file from the reference files workflow. The yaml option for this now takes `True` or `False` for whether to use genes instead of tiles
- Covariate files can their names to use in the model can now be given under the `covariates` section in the yaml. The name in the yaml will be used as the covariate's name in the model. If covariates are not to be used, the user can comment out the lines with the covariates names, but leave the overall covariates option so that it is `NULL`. I.e. just

```
covariates:
```

 with no lines indented after, or

```
covariates:
    #chrom_hmm: "/projects/dscott_prj/CCSRI_1500/exomes/ref/resources/E032_primaryBcells_fromPB_chromHMM_15state_segments.bed"
    #rep_timing: "/projects/dscott_prj/CCSRI_1500/exomes/ref/resources/RT_GM12878_Lymphocyte_Int90901931_hg19.bedgraph.gz"
```

## [1.2] - 2024-08-22

This update was authored by Sierra Gillis.

- Added an output directory for the prepare step
- Added more wildcards to the input maf symlink path for better knowing which was used for which run
- Was tested with an update made to `lcr-scripts/generate_smg_inputs` that uses threads to read in maf files faster, and only reads in relevant columns