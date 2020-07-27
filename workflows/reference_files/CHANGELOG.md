# Changelog

All notable changes to the `reference_files` subworkflow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.3] - 2020-07-17

This release was authored by Bruno Grande.

- Fixed bug with `download_dbsnp_vcf` where the output wasn't being redirected to the output file.

## [2.2] - 2020-07-16

This release was authored by Bruno Grande.

- The output reference files are no longer write-protected. 

## [2.1] - 2020-07-15

This release was authored by Bruno Grande.

- The log file for the STAR index creation rule was moved to outside of the STAR output directory, which gets deleted whenever the job fails.
- The provider for the Gencode GTF file was erroneously changed from "ucsc" to "ensembl". This version reverts that error. 

## [2.0] - 2020-06-29

This release was authored by Bruno Grande.

- This version added support for Sequenza reference files, namely the GC content wiggle (WIG) file and the dbSNP VCF file.
- The header of the `reference_files` Snakefile was stored in a separate file to increase user friendliness.
- Chromosomes Y and MT are excluded from the list of main chromosomes.

## [1.0] - 2020-05-02

This release was authored by Bruno Grande.

- This version can download genome FASTA files, Gencode annotations, and chromosome mappings between different providers (_e.g._ UCSC, Ensembl, NCBI).

- This version can index the genome FASTA file, create a BWA index, create a STAR index using the Gencode annotations for different overhang lengths, and create text and BED files for the main chromosomes.

- Chromosome names can automatically be converted between different providers using the [ChromosomeMappings GitHub repository](https://github.com/BrunoGrandePhD/ChromosomeMappings). For example, the Gencode annotations use the `chr` prefix, _i.e._ they use the UCSC chromosome names. This is stated in the Gencode download rule under `params` as `provider = "ucsc"`. For Ensembl genome builds, the chromosome names will automatically get converted from the UCSC style to the Ensembl style.

- All output files are converted to read-only as a safe guard once they're created by each rule.

- Downloads are hard linked (rather than symbolically linked) into the individual genome build directories to simplify permission management. When a hard link is created, it shares the same permissions as the target. In other words, if the target was made read-only, the hard link will automatically be read-only. Hard links should be avoided though when linking within the genome build directory to avoid issues where updating one rule would trigger other rules because the last-modified timestamp for all hard links gets updated if a change occurs in one of them. 

- The BWA and STAR indices capture the tool version in the file path in case of backwards-incompatible index changes. Because conda environments cannot be dynamic set (_e.g._ using functions), only one tool version can be used at a time. These are set in the `config/default.yaml`.

- The `create_bwa_index` rule doesn't actually symlink the FASTA file because that's technically not required. This way, the rule can avoid having to encode complicated logic to create a relative symlink since we can't use the `coreutils` conda environment to benefit from the `-r` (relative) option for `ln`. 
