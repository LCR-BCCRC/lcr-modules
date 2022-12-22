# Changelog

All notable changes to the `ecotyper` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-12-22

This release was authored by Kostia Dreval.

- The first version of this module only supports the bulk RNA-Seq for Lymphoma datasets. The support for carcinomas or single-cell data, as well as discovery of new ecotypes, can be implemented in later versions.
- The supplied gene expression matrix must be in wide format with first column containing gene names. The name of the first column is not checked and therefore can be anything, for example gene_symbol, Hugo_Symbol etc. By default, columns gene_id and ensembl_gene_id are dropped from matrix.
- The annotation data must contain column ID for sample ids. This is hard-coded in ecotyper. The annotation columns will be used to plot on the top of heatmaps. By default, column pathology is used. This can be specified in config. If no annotations to plot at heatmap, user shoul leave the config key empty (specify "").
- The ecotyper uses read.delim function to read the data, so unfortunatelly sample ids with special characters like - or . or ids starting with number are not properly handled. Therefore, this module implements pre- and post-processing rules to handle this.
- Ecotyper has some internal QC and therefore can filter out some samples during the run. Therefore, the state assignment label for each model may not always contain all the samples from input matrix. For more information about the sample filtering procedure please see the *Cell state quality control* section of the [EcoTyper paper] ([https://doi.org/10.1016/j.cell.2021.09.014](https://doi.org/10.1016/j.cell.2021.09.014)) methods
