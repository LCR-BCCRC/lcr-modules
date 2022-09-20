# Changelog

All notable changes to the `igcaller` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2022-05-10

This release was authored by Prasath Pararajalingam.

- Support paired/unpaired mode
- IgCaller only available on github, therefore require user to clone IgCaller repo locally and set path to the repo in config.yaml
- IgCaller repo contains required reference files separated by human genome reference versions. -V parameter of IgCaller takes either 'hg19' or 'hg38' to choose which reference version to use. Parameterized possible values for this parameter in config.yaml under "genome_version" and key-value pairs corresponding to various genomes. Up to the user to add any new genomes.
- Similar to above, -C takes either 'ucsc' or 'ensembl' to decide whether to use or omit 'chr'-prefix. Parameterized this option in config.yaml under "chr_annotation" with key-value pairs corresponding to various genomes.
- -seq parameter chooses sequencing type (i.e., wgs, wes, or capture). Used switch on seq_type wildcard to choose either wgs for genomes or wes for capture. 
