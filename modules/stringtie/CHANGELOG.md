# Changelog

All notable changes to the `stringtie` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-11-30

This release was authored by Krysta Coyle.

- awk script necessary to introduce XS tags into STAR-aligned BAMs. Without XS tags, stringtie generates single-exon transcripts.
- Reference GTF mandatory, highly recommended for well-annotated genomes.
- Initial estimates of memory are conservative.
- Module designed to work with STAR bams output from modules/STAR/1.4.
