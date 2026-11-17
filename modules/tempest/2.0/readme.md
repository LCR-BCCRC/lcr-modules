<p align="center"><img height="300" alt="TEMPEST_logo" src="https://github.com/user-attachments/assets/ff11d9d8-1e14-40dc-a6c8-b7a2f557ce1a" />

# TEMPEST v2.0

TEMPEST (Temporal Evaluation of Mutations in Plasma cfDNA with Error-suppressed Somatic Tracking) is a pipeline designed
to process hybrid capture sequencing samples from liquid biopsies of cancer patients. The pipeline can account
for serial sampling of the same patient and track variants across time points. The pipeline relies on a duplex consensus deduplication
step using UMIs. Version 2 of TEMPEST follows fgbio's suggested workflow for duplex consensus deduplication, and is more suitable for large
input files.

The pipeline consists of two main workflows: one aligning FASTQs and deduplicating reads using UMIs,
and then variant calling and filtering using SAGE (https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md) as the caller.

There is an optional workflow to compile per-patient summary reports. The example code given here is a rough template, and will likely require
the creation of a custom Jupyter notebook template for your project.

The pipeline requires Conda and Snakemake, and was built using Snakemake v8.10.8.

## Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Read Support Filters](#read-support-filters)
    - [Hotspot and Notspot Lists](#hotspot-and-notspot-lists)
    - [Phase Sets](#phase-sets)
    - [UMI Support](#umi-support)
    - [augmentMAF](#augmentmaf)
- [Artifact Filtering](#artifact-filtering)
- [Resources](#resources)
- [Samplesheet requirements](#samplesheet-requirements)
- [Example Usage](#example-project-snakemake-workflow)
- [Config parameters](#config-parameters)
- [Patient Reports](#note-about-patient-reports)
- [Panel of Normals](#panel-of-normals)
- [DLBCL Panel Hotspot and Notspot Lists](#dlbcl-panel-hotspot-and-notspot-lists)

# Pipeline Overview
TEMPEST is designed to take raw FASTQ files and produce QC and variant call files. The primary workflow for variant calling is designed
to utilize matched normal samples from a patient, using the workflow `.../variantTEMPEST.smk`, and filtering has been optimized for
using this workflow. However, unpaired samples can have variants called by using `.../unpaired_variantTEMPEST.smk`. However, besides
using gnomAD frequencies and a PON for filtering, no other steps attempt to remove germline mutations.

<p align="center"><img height="800" alt="TEMPEST pipeline (1)" src="https://github.com/user-attachments/assets/2bfca34f-6e51-4ba4-b3a1-3b7f6a3fe73e" />

## Read Support Filters
Final variant filtering is carried out by `.../utils/custom_filters.py`. To help determine why variants were included in the final output, a
`variant_source` column is added to the final .maf, which lists the reasons a variant passed filtering, possible values are: high_vaf, hotspot,
additional_maf or phase_group. 

It should be noted that part of TEMPEST's approach to filtering variants is built on the premise that there can be some contamination of tumour
DNA in the matched normal. Hypothetically, this could occur as a result of poor fractionation of a plasma sample or from circulating tumour cells
being included with PBMCs. Regardless of the mechanism, when developing TEMPEST on lymphoma patient samples it was repeatedly found that read support 
can be found in matched normal samples for tumour variants, including known hotspot variants. Therefore, some flexibility for alt reads in the matched
normal was built into TEMPEST. This can be adjusted in the config.

Values shown are defaults which can be specified in the config. <br>
<p align="center"><img height="600" alt="Variant Filtering (2)" src="https://github.com/user-attachments/assets/b1b501f6-3baa-404f-80f3-3ee9b5f62a77" />


### augmentMAF

TEMPEST uses a script called `.../utils/augmentMAF.py` to augment final variant calls in a given sample with still detectable, but uncalled, variants found in other samples from the same patient. It does this by searching for read support for variants from all available samples from a patient, in the given index sample .bam file. These augmented variants provide increased sensitivity, particularly when tracking variants temporally.

As variants being augmented have already passed filtering in their origin sample, criteria for their augmentation is less strict. Included in this augmentation process, is the tracking of phased variants. The augmentMAF script combines intersecting phased variant sets from all input samples (index and additional), then determines if these variants are phased in the given sample. Variants that were found to be phased are given a value of TRUE in the column `is_phased` added to the final output .maf.

For augmented variants, their originating sample(s) are listed in the `origin_samples` column added to the final output .mafs.

### Hotspot and Notspot Lists

TEMPEST uses custom hotspot and notspot (aka blacklist) positions in conjunction with other read support parameters in a custom script. Therefore, user-provided lists are integral to enhancing performance of filtering. The lists should be provided as a txt file with no header, with each line being a single position, formatted like: `chr:position`.

### Phase Sets

Phasing information is used in two places in the pipeline: in the the read support filtering and when augmenting the final variant set (described above).

When calling variants SAGE assigns variants with overlapping read evidence to local phase sets, and this information is preserved and used
by TEMPEST when doing its own read support filtering (diagram above). One or more variants from a phase set need to be present for TEMPEST to consider
variants phased.

### UMI Support
Using a custom script `.../utils/FetchVariantUMIs.py`, TEMPEST annotates variant calls in each .maf file with metrics related to the UMI family
sizes of reads containing alt alleles.

- **UMI_mean**: the mean UMI family size for all alt allele reads
- **UMI_max**: the maximum UMI family size amongst all alt allele reads
- **UMI_3_count**: How many alt allele reads have a UMI family size of ≥ 3 (ie had some error correction)

For filtering UMI_3_count is used to add requirements for error correction on reads that support alt alleles.

### CHIP Variants

TEMPEST flags, but does not remove, variants that could be a result of clonal hematopoiesis of intermediate potential (CHIP) with a boolean column `CHIP` in the final .maf.

Variants are flagged if: 1) the alt allele VAF is ≥ 2% in the matched normal or 2) the alt allele vaf is x3 larger in the matched normal than in the tumour sample. 

However, TEMPEST is more strict with variants in commonly reported genes that accumulate age-related CHIP mutations: "DNMT3A","TET2","ASXL1","PPM1D","TP53","JAK2","SF3B1","SRSF2", and flags any variant in these genes with an alt_read_count > 0 in the matched normal.

# Artifact Filtering

Technical artifacts, in the form of low VAF variant calls, are an inherent part of NGS assays. TEMPEST takes advantage of UMI error correction to surpress this noise, but further methods are also availble. 

1) Custom Notspot/Blacklists: user provided known trouble spots to be filtered out
2) Panel of Normals: a user created panel of sites commonly mutated in normal samples sequenced using the same assay intended to be used on tumour/experimental samples.
3) Background mutation rates: Calculate the background mutation rate and SD for all positions in the target space of a hybrid capture panel, and require the VAF of a candidate variant to be a specified times greater than the mean (default 2xSD).

Assistance in how to make a PON and a background mutation rate index are detailed below.

## Backgroud Mutation Rate Index

TEMPEST can take an input mutation rate index (.tsv file) and filter variants based on how noisy a position is. The background mutation rate can be calculated by using a selection of input samples sequenced with your assay, ideally normal or healthy samples. To can use the `artifact_alert` lcr-module to create this index, or make your own solution. 

Either way here is en example of the required format:

chromosome | position | ref_base | n_samples | A_mean    | A_std      | T_mean    | T_std      | C_mean    | C_std      | G_mean    | G_std      | INS_mean | INS_std    | DEL_mean | DEL_std    |
|------------|----------|----------|-----------|-----------|------------|-----------|------------|-----------|------------|-----------|------------|----------|------------|----------|------------|
| chr1       | 2556665  | A        | 358       | 0.0       | 0.0        | 1.031e-05 | 0.00019502 | 4.773e-05 | 0.00038169 | 3.66e-05  | 0.00029125 | 0.0      | 0.0        | 0.0      | 0.0        |
| chr1       | 2556666  | T        | 358       | 5.05e-06  | 9.557e-05  | 0.0       | 0.0        | 0.00023905| 0.00290158 | 9.987e-05 | 0.0010797  | 6.15e-06 | 0.00011641 | 0.0      | 0.0        |
| chr1       | 2556667  | G        | 358       | 0.00011274| 0.00066822 | 6.775e-05 | 0.00051612 | 6.379e-05 | 0.00088875 | 0.0       | 0.0        | 0.0      | 0.0        | 5.009e-05| 0.00070783

## Panel of Normals

To help exclude technical artifacts unique to a hybrid capture panel, it is recommended that you sequence some cfDNA
from "healthy" individuals and use their data to create a panel of normals. That VCF is required by the pipeline for filtering,
but if you don't have one you can provide an empty dummy file.

To assist in the creation of a PON, a workflow using GATK is included in this module: `make_somaticPON.smk`

<p align="center"><img height="400" alt="TEMPEST pipeline" src="https://github.com/user-attachments/assets/f6877543-ad99-416d-ad42-59526c92b730" />

# Resources

## SAGE
Variant calling is done initially using SAGE:
https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md

If you need to download resources for SAGE, check their downloads page:
https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md#variants

# Samplesheet requirements
The pipeline requires a user-created samplesheet with details about the samples.

| Field           | Allowed values/format                    | Description/notes                                                                                 |
|-----------------|------------------------------------------|----------------------------------------------------------------------------------------------------|
| sample_id       | string                                   | Must be unique                                                                                    |
| patient_id      | string                                   | Must match across all samples from the same patient                                               |
| tissue_status   | tumor, normal                            | Sample tissue classification                                                                      |
| matched_normal  | sample_id, or "unmatched"                | The sample_id of the matched normal sample, or "unmatched" if using the unmatched workflow       |
| genome_build    | hg37, hg38                               | Ensure the corresponding resource versions are provided in the config                             |
| seq_type        | capture, genome                          | Sequencing method used                                                                             |
| capture_space   | string (panel name)                      | Name of the capture panel; must correspond to a dictionary entry in the config                     |

# Example project Snakemake workflow

Here is a brief example of how to set up your project Snakemake workflow and include the necessary inputs. This
assumes you have also completed the config file.

someworkflow.smk
```
import pandas as pd

# set up a config and samplesheet, ensure all required columns are present
configfile: path/to/completed/config.yaml
SAMPLESHEET = pd.read_csv(samplesheet.tsv, sep="\t")

# add sample table into config
config["lcr-modules"]["_shared"]["samples"] = SAMPLESHEET.copy()

# include workflows you want to use
include: .../lcr-modules/modules/tempest/1.0/align_and_deduplicateTEMPEST.smk
include: .../lcr-modules/modules/tempest/1.0/variantTEMPEST.smk

# specify your rule all
rule all:
    input:
        rules.all_bams.input,
        rules.all_sage.input,
    threads: 1
    resources:
        mem_mb = 500
```

The workflows include functions for dynamically setting memory requirements when submitting jobs to a compute cluster.

# Config parameters

The config file contains many values, with defaults for variant calling and filtering. Those defaults were optimized
for cfDNA samples from DLBCL patients, with mean target coverage between 400–1000. 

Some key config options for variant calling and filtering are detailed here:

Config Option                | Description                                                                                                  | Default |
|-----------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| mask_threshold              | Minimum fraction of uncalled (N) bases at a position to remove it                                            | 0.05    |
| mask_count                  | Minimum number of N bases required to mask a position                                                         | 8       |
| hard_min_vaf                | Hard minimum VAF floor; should be ≤ SAGE's tier-specific cutoffs to speed processing                          | 0.002   |
| tumor_panel_min_vaf         | First-pass VAF threshold given to SAGE for panel positions                                                     | 0.005   |
| novel_vaf                   | VAF threshold for novel variants; only hotspots and phased variants may be lower                            | 0.01    |
| panel_max_germ_rel_raw_bq   | Max proportion of good-quality alt reads allowed in normal vs tumor (cumulative BQ); SAGE arg (SAGE default 0.04) | 0.08    |
| hotspot_max_germ_rel_raw_bq | Max proportion for hotspot positions (cumulative BQ); SAGE arg (SAGE default 0.25)                          | 0.08    |
| hotspot_max_germline_vaf    | Max VAF allowed for germline variants in hotspots                                                           | 0.1      |
| min_germline_depth          | Minimum depth in the normal (germline) sample                                                               | 50      |
| max_normal_alt_depth        | Maximum alt depth allowed in the normal; hard filter not applied at hotspot positions                       | 10      |
| min_map_qual                | Minimum mapping quality                                                                                    | 40      |
| min_alt_depth               | Minimum alt read count                                                                                      | 5       |
| min_t_depth                 | Minimum tumor depth                                                                                          | 100     |
| exac_max_freq               | Maximum population frequency in gnomAD to retain a variant                                                    | 0.01    |
| min_UMI_3_count             | Minimum UMI_3 count for variants with t_alt_count > low_alt_thresh                                          | 2       |
| low_alt_thresh              | Alt-count threshold below which stricter UMI_3 rules apply                                                  | 50      |
| low_alt_min_UMI_3_count     | Minimum UMI_3 count for low-alt-count variants (t_alt_count ≤ low_alt_thresh)                               | 3       |
| aug_min_UMI_3_count         | Minimum UMI_3 count used in augmented MAF filter criteria                                                   | 1       |
| phase_id_col                | SAGE phase/ID column name                                                                                   | "LPS"   |
| phased_min_t_alt            | Minimum tumor alt count for phased variants                                                                 | 3       |
| background_rates            | Path to file with background mutation rate index                                                            | None     |
| min_background_samples      | Min number of samples used to calculate background mutation rate to allow a position to be passable          | 20     |
| background_n_std            | number of SD above the mean the VAF of a candidate variant needs to be to pass filter                        | 2       |
| max_background_rate         | hardcap on background mutaion rate. All positions with this rate or higher are filtered out!                 | 0.1     |


# Note about patient reports

Due to the uniqueness of the required analyses of different projects `patient_reports.smk` is not a one-click solution
for all projects. Therefore, the workflow and code here will likely need to be adapted for your usage. Also included here is
an example templated Jupyter notebook that will be used to compile a report `report_template.ipynb`, and the accompanying 
python script `compile_report.py` that completes the report creation.

These scripts can be found:
```
.../lcr-modules/modules/Tempest/1.0/patient_reports/report_template.ipynb
.../lcr-modules/modules/Tempest/1.0/patient_reports.smk
.../lcr-modules/modules/Tempest/1.0/patient_reports/compile_report.py
```

The required config entries for this example workflow are:

```
cfDNA_patient_reports:
    report_template: ".../ctDNA_Pro_Pipeline/patient_reports/report_template.ipynb"
    include_lymphgen: true # boolean

```
As this code was developed on lymphoma patient samples, it can integrate the LymphGen classifier outputs as implemented in`lcr-modules`. You can simply
set that to `false` in the config if you do not need it.

# DLBCL Panel Hotspot and Notspot Lists

This pipeline was designed using the [DLBCL-STORM](https://github.com/morinlab/DLBCL-STORM "DLBCL-STORM GitHub repo") hybrid capture panel,
which has its own custom hotspot/notspot lists.
