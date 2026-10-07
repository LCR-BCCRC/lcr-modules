# Changelog

All notable changes to the `nanovote` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2026-09-14

This release was authored by Giuliano Banco.

NanoVote is an ensemble somatic variant-calling module that combines calls from ClairS, ClairS-TO, and DeepSomatic into a consensus call set.

Orchestrates the three caller submodules and merges their outputs.

Uses Starfish to identify variants supported by multiple callers and produces consensus, union, and 2+ output VCFs.

The 2+ callset contains all variants supported by at least two callers. This set is less strict when calling variants and is likely to produce more false positives than the final set.

The final callset contains all variants supported by ClairS-TO and at least one other caller. This set is more strict when calling variants and reduces false positives at the cost of producing more false negatives than the 2+ set.
