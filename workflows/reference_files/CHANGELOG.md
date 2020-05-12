# Changelog

All notable changes to the `reference_files` subworkflow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-05-02

This release was authored by Bruno Grande.

- This version can download genome FASTA files, Gencode annotations, and chromosome mappings between different providers (_e.g._ UCSC, Ensembl, NCBI).

- This version can index the genome FASTA file, create a BWA index, create a STAR index using the Gencode annotations for different overhang lengths, and create text and BED files for the main chromosomes.

- Chromosome names can automatically be converted between different providers using the [ChromosomeMappings GitHub repository](https://github.com/BrunoGrandePhD/ChromosomeMappings). For example, the Gencode annotations use the `chr` prefix, _i.e._ they use the UCSC chromosome names. This is stated in the Gencode download rule under `params` as `provider = "ucsc"`. For Ensembl genome builds, the chromosome names will automatically get converted from the UCSC style to the Ensembl style.

- All output files are converted to read-only as a safe guard once they're created by each rule.

- Everything is hard linked (rather than symbolically linked) to simplify permission management. When a hard link is created, it shares the same permissions as the target. In other words, if the target was made read-only, the hard link will automatically be read-only.

- The BWA and STAR indices capture the tool version in the file path in case of backwards-incompatible index changes. Because conda environments cannot be dynamic set (_e.g._ using functions), only one tool version can be used at a time. These are set in the `config/default.yaml`.
