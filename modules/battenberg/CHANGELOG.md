# Changelog

All notable changes to the `battenberg` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-15

This release was authored by Ryan Morin.

I modified the Battenberg source code to take an optional "chr_prefixed_genome" argument and some modifications to enable battenberg to work with genomes that have "chr" prefixes in the chromosome names. The installation uses conda to create an R installation with all required prerequisites with the exception of ASCAT and Battenberg, which are installed from github. 
This installation is currently accomplished by a _install_battenberg rule but there is probably a better way (subworkflow?).  



- No module design decisions explained here yet.
