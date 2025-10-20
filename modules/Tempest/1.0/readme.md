# TEMPEST

TEMPEST (Temporal Evaluation of Mutaions in Plasma cfDNA with Error-suppressed Somatic Tracking) is a pipeline designed
to process hybrid capture sequencing samples from liquid biopsies of cancer patients. The pipeline is able to take into account
serial sampling of the same patient, and track variants across timepoints. The pipeline relies on a duplex consensus deduplication
step using UMIs.

The pipeline consists of two main workflows, one aligning fastqs and deduplicating reads using UMIs,
and then variant calling and filtering using SAGE [https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md] as a caller.

There is an optional workflow to compile per-patient summary reports. The example code given here is a rough template, and will likely require
the creation of a custom jupyter notebook template for your project.

The pipeline uses snakemake (v. 8.10.8) and conda.

# Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Resources](#resources)
- [Samplesheet requirements](#samplesheet-requirements)
- [Example Usage](#example-project-snakemake-workflow)
- [Config parameters](#config-parameters)
- [Patient Reports](#note-about-patient-reports)
- [Panel of Normals](#panel-of-normals)
- [DLBCL Panel, hotstpots sand notspots](#dlbcl-panel-hotspot-and-notspot-lists)

# Pipeline Overview


# Resources 

## SAGE
Variant calling is done initially using SAGE:
https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md

If you need to download resources for SAGE check their downloads page:
https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md#variants

# Samplesheet requirements
The pipeline requires a user created samplesheet with details about the samples.

| Field           | Allowed values/format                    | Description/notes                                                                                 |
|-----------------|------------------------------------------|----------------------------------------------------------------------------------------------------|
| sample_id       | string                                   | Must be unique.                                                                                    |
| patient_id      | string                                   | Must match across all samples from the same patient.                                               |
| tissue_status   | tumor, normal                            | Sample tissue classification.                                                                      |
| matched_normal  | sample_id, or "unmatched"                | The sample_id of the matched normal sample, or "unmatched" if using the unmatched workflow.       |
| genome_build    | hg37, hg38                               | Ensure the corresponding resource versions are provided in the config.                             |
| seq_type        | capture, genome                          | Sequencing type (use genome for lymphgen).                                                         |
| capture_space   | string (panel name)                      | Name of the capture panel; must correspond to a dictionary entry in the config.

# Example project snakemake workflow

Here is a brief example of how to setup your project snakemake workflow and include the necessary inputs. This
assumes you have also completed the config file.

someworkflow.smk
```
import pandas as pd

# setup config and samplesheet
configfile: path/to/completed/config.yaml
SAMPLESHEET = pd.read_csv(samplesheet.tsv, sep="\t")
# add sample table into config
config["lcr-modules"]["_shared"]["samples"] =  SAMPLESHEET.copy()
############################ include any required filtering of samples to be processed, or add column checks

# include workflows you want to use
include: .../lcr-modules/modules/Tempest/1.0/cfdna_UMI_workflow.smk
inlcude: .../lcr-modules/modules/Tempest/1.0/cfDNA_SAGE_workflow.smk

# specify your rule all

rule all:
    input:
        rules.all_bams.input,
        rules.all_sage.input,
    threads: 1
    resources:
        mem_mb = 500

```

# Config parameters

The config file contains many values, with defaults for variant calling and filtering. Those defaults were optimized
for cfDNA samples from DLBCL patients, with mean target coverage between 400-1000. 

Some key config options for variant calling and filtering are detailed here:

Config Option                | Description                                                                                                  | Default |
|-----------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| mask_threshold              | Minimum fraction of uncalled (N) bases at a position to remove it                                            | 0.05    |
| mask_count                  | Minimum number of N bases required to mask a position                                                         | 8       |
| hard_min_vaf                | Hard minimum VAF floor; should be ≤ tier-specific cutoffs to speed processing.                               | 0.002   |
| tumor_panel_min_vaf         | First-pass VAF filter on panel positions.                                                                    | 0.005   |
| novel_vaf                   | VAF threshold for novel variants; only hotspots and phased variants may be lower.                            | 0.01    |
| panel_max_germ_rel_raw_bq   | Max proportion of good-quality alt reads allowed in normal vs tumor (cumulative BQ); SAGE arg (SAGE default 0.04). | 0.08    |
| hotspot_max_germ_rel_raw_bq | Max proportion for hotspot positions (cumulative BQ); SAGE arg (SAGE default 0.25).                          | 0.08    |
| min_germline_depth          | Minimum depth in the normal (germline) sample.                                                               | 50      |
| max_normal_alt_depth        | Maximum alt depth allowed in the normal; hard filter not applied at hotspot positions.                       | 10      |
| min_map_qual                | Minimum mapping quality.                                                                                     | 40      |
| min_alt_depth               | Minimum alt read count (post-filtering of index variants).                                                   | 5       |
| min_t_depth                 | Minimum tumor depth (especially important for unpaired workflows).                                           | 100     |
| exac_max_freq               | Maximum population frequency in ExAC to retain a variant.                                                    | 0.01    |
| min_UMI_3_count             | Minimum UMI_3 count for variants with t_alt_count > low_alt_thresh.                                          | 2       |
| low_alt_thresh              | Alt-count threshold below which stricter UMI_3 rules apply.                                                  | 50      |
| low_alt_min_UMI_3_count     | Minimum UMI_3 count for low-alt-count variants (t_alt_count ≤ low_alt_thresh).                               | 3       |
| aug_min_UMI_3_count         | Minimum UMI_3 count used in augmented MAF filter criteria.                                                   | 1       |
| phase_id_col                | SAGE phase/ID column name.                                                                                   | "LPS"   |
| phased_min_t_alt            | Minimum tumor alt count for phased variants.                                                                 | 3       |


# Note about patient reports
Due to the uniquness of the required analyses of different projects the patient_reports.smk is not a one click solution
for all projects. Therefore the workflow and code here will likely need to be adapted for your usage. 

To run it you need to provide in the config (under lcr-modules) the following:

```
cfDNA_patient_reports:
    report_template: ".../ctDNA_Pro_Pipeline/patient_reports/report_template.ipynb"
    include_lymphgen: True # boolean

```
As this code was developed on lymphoma patient samples, it has the option to integrate the lymphgen classfier lcr-module outputs. You can simply
set that to "False" in the config if you do not need it.


Which includes the path to a jupyter notebook template. To be able to full understand what is needed to get this workflow going
for you, you need to look at both the example template, the reports snakemake, and the script that compiles the reports:

.../lcr-modules/modules/Tempest/1.0/patient_reports/report_template.ipynb
.../lcr-modules/modules/Tempest/1.0/patient_reports.smk
.../lcr-modules/modules/Tempest/1.0/patient_reports/compile_report.py

# Panel of Normals

To help exclude technical artifacts unique to a hybrid capture panel, it is recommended you sequence some cfDNA
from "healthy" individuals, and use their data to create a panel of normals. That VCF is required by the pipeline for filtering,
but if you dont have one you can give it a empty dumby file.

To assist in the creation of a PON a workflow using GATK is included in this module: make_somaticPON.smk.

# DLBCL Panel Hotspot and Notspot lists

This pipeline was designed during a project using [DLBCL-STORM](https://github.com/morinlab/DLBCL-STORM "DLBCL-STORM GitHub repo"),
which has its own custom hybrid capture panel and hotspot/notspot lists.