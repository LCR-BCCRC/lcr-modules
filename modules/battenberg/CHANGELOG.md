# Changelog

All notable changes to the `battenberg` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-17

This release was authored by Ryan Morin.

I modified the Battenberg source code to take an optional "chr_prefixed_genome" argument and some modifications to enable battenberg to work with genomes that have "chr" prefixes in the chromosome names. The installation uses conda to create an R installation with all required prerequisites with the exception of ASCAT and Battenberg, which are installed from github. 
This installation is currently accomplished by a _install_battenberg rule but there is probably a better way (subworkflow?).  

In contrast to the Sequenza module, I haven't embedded cnv2igv.py in the module. Instead, the user must specify the path to their LCR-scripts repository. Similar to my original version of Sequenza, I relied on the calc_sex_status.sh script because it was working. More elegant solutions are welcome. 
