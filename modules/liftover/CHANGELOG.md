# Changelog

All notable changes to the `liftover` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0] - 2021-12-29

This release was updated by Laura Hilton

- Add support for bedpe files

## [1.2] - 2021-01-11

- This version significantly improves module functionality by introduction of the following changes: support of biderectional conversion of genomic coordinates (hg19 -> hg38 and hg38 -> hg19), sorting of converted bed files, filtering converted outputs to include standard chromosomes only, and filling of the gaps with empty segments to ensure compatiability with Level 3 analyses. In addition, chain files is now generated through reference_files workflow and therefore there is no need to include it with the module `src`.


## [1.1] - 2020-09-05

- This version includes 2 important updates: now the module can filter samples automaticaly to lift only 
  samples in hg38 genome build, which makes module more general and involves less input from user. In addition,
  the ability to add parent directory is now avaliable in config file, which enables tracking of not only the 
  variant caller, but also version of the module used.


## [1.0] - 2020-07-01

This release was authored by Kostiantyn Dreval.

- Initial draft of liftover module. It uses custom script that does automatic conversion of input file 
  from seg to bed format and vice versa. To preserve the name of the variant caller that produced intial 
  seg file, the variant caller information is handled as `{tool}` wildcard throughout the module and is 
  parsed to the name of lifted files. Since bed file, unlike files in seg format, does not contain header,
  the header of each input file supplied for conversion by liftover is stored as a separate file, which can
  be made temporary in later versions of te module.
