# Changelog

All notable changes to the `fusioncatcher` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2025-01-17

This release was authored by Ryan Morin.

* The [fusioncatcher](https://github.com/ndaniel/fusioncatcher/tree/master) pipeline conveniently has a conda recipe available, which allowed almost all the dependencies to be handled automatically. 
* I haven't incorporated the installation of the reference files into this initial version of the module. Instead, the module expects the user to do that themselves. Because the process is straightforward, well documented and doesn't require any dependencies, this seemed reasonable to me.
* The conda environment has one reproducible defect. The faToTwoBit binary complains about a missing .so file. The solution was to re-download that file in a post-deploy script and replace the broken file. 
