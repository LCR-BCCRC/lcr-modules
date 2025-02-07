# cfDNA Pipeline documentation

This pipeline coordinates the analysis of DNA sequencing data from hybrid capture assays of cfDNA. 

The pipeline consists of two main workflows, one aligning fastqs and deduplicating reads using UMIs,
and then variant calling and filtering using SAGE as a caller.

The pipeline requires some resources not contained in this repo, links for some are provided below.

For GSC users these are probably floating around already.

The pipeline uses snakemake and conda.

# Resources 

## SAGE
https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md

If you need to download resources for SAGE try downloading the prepackaged set of them that they provide from:
https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline/


# Samplesheet requirements
The samplesheet has the following required columns

sample_id: must be unique
patient_id: must be the same as other samples from the same patient
tissue_status: tumor or normal
matched_normal: the sample_id of the matched normal sample. OR a value of "unpaired" if using the unmatched sample workflow.