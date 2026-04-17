# Changelog

All notable changes to the `sniffles` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1] - 2026-04-17

This release was authored by Giuliano Banco.

In this release, the sniffles conda environment was updated to use sniffles version 2.7.5.
Sniffles flags can now be supplied in the config, establishing support for custom flag use.
The deprecated flag --non-germline was removed.
To support the use of cram files, the input rule now symlinks the index file as a crai.
