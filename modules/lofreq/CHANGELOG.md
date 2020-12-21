# Changelog

All notable changes to the `lofreq` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-12-21

- add .crai index for cram support

## [1.0] - 2020-12-19

- Previous bug fix didn't handle some of the situations in which the non-ACTG reference bases causing garbled VCFs. This bundled lofreq binary avoids this situation entirely.
- This will probably cause issues with compatability on other platforms but this is better than our current version that fails consistently on some samples. 

## [1.0] - 2020-12-04

This release was authored by Ryan Morin 

This fixes a bug that was causing lofreq to die consistently on some cases. In the process of setting up the snakefile to use the bundled version of the patched script, I changed the structure of the src directory. This simplifies how the pipeline calls other bundled scripts as well. I'm considering this approach as a model for how to do this in other modules that use bundled helper scripts. 

## [1.0] - 2020-11-02

This release is authorized by Kostia Dreval

- It was found that if lofreq run is killed by any reasom midway through and restarted to continue on same samples, it is exiting with error refusing to overwrite the files of partially completed outputs. To prevent this from happening, this release features added conditional statement before launching the lofreq, that checks whether target outputs are present.
- The filtering script was modified to replace the `AF` field in output vcf files with `VAF` to ensure consistency with GNOMAD.


## [1.0] - 2020-10-27

This release is authorized by Kostia Dreval

- It has been found that in the lofreq outputs are present non-stabdard chromosomes and non-ATCG characters. In this modification is introduced
  bash script associated with the module that will filter out these features. Because the filtered vcf file needs to be indexed, the output of the script should be bgziped for future use with bcftools tabix. In addition, the naming scheme of final outputs is modified to be consistent with other
  variant callers.
- At current version lofreq is not calling indels, as it requires the input bam files to be pre-processed. The option to call indels will be
  introduced in later versions.

## [1.0] - 2020-10-15

This is a release of updated initial draft from Bruno on 2020-06-30.

- Kostia have removed the `op.retry` function for the `mem_mb` resource in `_lofreq_run` rule. Also, added `bam` option in resources to allow for 
  restriction of bam files processed simultaneously. Increased number of threads and memory resources comparing to initial draft to make sure 
  genomes are running on a cluster. Tested the module on both genomes and exome sample from demo workflow.
