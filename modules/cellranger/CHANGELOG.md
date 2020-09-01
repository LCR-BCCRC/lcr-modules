# Changelog

All notable changes to the `Cell Ranger` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0] - 2020-04-25

This module is specialized for single cell mRNA sequencing analysis

The `samples.tsv` file for this analysis must include the following columns: `sample_id`, `index`, `lane`, `chip_id`, `analysis` and `seq_type`.
  - the columns correspond to the required columns in Cell Ranger samplesheet format
  - `analysis` column should be filled by either `count` or `vdj` which corresponds to which Cell Ranger analysis is run on it after    `Cell Ranger count`
  - `seq_type` must be specified as sc_mrna

Since the Cell Ranger software is self-contained, conda environments are not utilized to run it. To use a version of Cell Ranger, users can modify `CFG["software"]` or `config["lcr-modules"]["Cell Ranger"]["software"]`, to contain a path to the Cell Ranger executable

The default software is the most recent version, Cell Ranger v4.0.0

Stamps are used as dummy inputs and outputs file to determine rule dependencies.


