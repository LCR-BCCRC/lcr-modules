# Changelog

All notable changes to the `slms_3` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-01-11

This release was authored by Laura Hilton.

- SLMS-3: SNVs and indels called by three or more of Strelka2, LoFreq, Mutect2, and SAGE. 
- This module incorporates all the variant calling steps for the SLMS-3 pipeline, which can be used for FFPE and FF tumours run with a matched or unmatched normal. 
- The versions of each submodule can be specified in the config. 
- If running this pipeline on a large number of samples, resource throttling is available. Currently each high resource use job is assigned the following resource values: 
    ``` 
    strelka-1.1: bam=1
    manta-2.3: bam=1
    sage: sage=1
    lofreq: lofreq=1
    mutect2: mutect_chrom=1, mutect_pup=1
    ````

    Therefore you might consider launching Snakemake with `--resources bam=40 sage=20 lofreq=20 mutect_chrom=100 mutect_pup=20` which will ensure a max of 40 Strelka/Manta jobs running simultaneously, max 20 SAGE jobs, etc. If there is a large number of unmatched normal samples in your pipeline, you may need to restrict this further to prevent I/O bottlenecks. 

    
