# TEMPEST

TEMPEST (Temporal Evaluation of Mutations in Plasma cfDNA with Error-suppressed Somatic Tracking) is a pipeline designed
to process hybrid capture sequencing samples from liquid biopsies of cancer patients. The pipeline can account
for serial sampling of the same patient and track variants across time points. The pipeline relies on a duplex consensus deduplication
step using UMIs.

The pipeline consists of two main workflows: one aligning FASTQs and deduplicating reads using UMIs,
and then variant calling and filtering using SAGE (https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md) as a caller.

There is an optional workflow to compile per-patient summary reports. The example code given here is a rough template, and will likely require
the creation of a custom Jupyter notebook template for your project.

The pipeline uses Snakemake (v8.10.8) and Conda.

# Table of Contents
- [Pipeline Overview](#pipeline-overview)
- [Variant Filtering Overview](#variant-filtering-overview)
- [Resources](#resources)
- [Samplesheet requirements](#samplesheet-requirements)
- [Example Usage](#example-project-snakemake-workflow)
- [Config parameters](#config-parameters)
- [Patient Reports](#note-about-patient-reports)
- [Panel of Normals](#panel-of-normals)
- [DLBCL Panel Hotspot and Notspot Lists](#dlbcl-panel-hotspot-and-notspot-lists)

# Pipeline Overview

<img alt="TEMPEST pipeline" src="https://github.com/user-attachments/assets/ee7a4e1d-88e1-451f-ad79-c3849d350a2e" height="800" />

## Variant Filtering Overview

<img alt="Variant Filtering" src="https://github.com/user-attachments/assets/046d1def-6052-4eaf-92ab-365b62b8ecd9" height="600" />

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
include: .../lcr-modules/modules/Tempest/1.0/cfdna_UMI_workflow.smk
include: .../lcr-modules/modules/Tempest/1.0/cfDNA_SAGE_workflow.smk

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
for cfDNA samples from DLBCL patients, with mean target coverage between 400–1000. 

Some key config options for variant calling and filtering are detailed here:

Config Option                | Description                                                                                                  | Default |
|-----------------------------|--------------------------------------------------------------------------------------------------------------|---------|
| mask_threshold              | Minimum fraction of uncalled (N) bases at a position to remove it                                            | 0.05    |
| mask_count                  | Minimum number of N bases required to mask a position                                                         | 8       |
| hard_min_vaf                | Hard minimum VAF floor; should be ≤ SAGE's tier-specific cutoffs to speed processing                          | 0.002   |
| tumor_panel_min_vaf         | First-pass VAF threshold given to SAFE for panel positions                                                     | 0.005   |
| novel_vaf                   | VAF threshold for novel variants; only hotspots and phased variants may be lower                            | 0.01    |
| panel_max_germ_rel_raw_bq   | Max proportion of good-quality alt reads allowed in normal vs tumor (cumulative BQ); SAGE arg (SAGE default 0.04) | 0.08    |
| hotspot_max_germ_rel_raw_bq | Max proportion for hotspot positions (cumulative BQ); SAGE arg (SAGE default 0.25)                          | 0.08    |
| min_germline_depth          | Minimum depth in the normal (germline) sample                                                               | 50      |
| max_normal_alt_depth        | Maximum alt depth allowed in the normal; hard filter not applied at hotspot positions                       | 10      |
| min_map_qual                | Minimum mapping quality                                                                                    | 40      |
| min_alt_depth               | Minimum alt read count                                                                                      | 5       |
| min_t_depth                 | Minimum tumor depth                                                                                          | 100     |
| exac_max_freq               | Maximum population frequency in gnomaD to retain a variant                                                    | 0.01    |
| min_UMI_3_count             | Minimum UMI_3 count for variants with t_alt_count > low_alt_thresh                                          | 2       |
| low_alt_thresh              | Alt-count threshold below which stricter UMI_3 rules apply                                                  | 50      |
| low_alt_min_UMI_3_count     | Minimum UMI_3 count for low-alt-count variants (t_alt_count ≤ low_alt_thresh)                               | 3       |
| aug_min_UMI_3_count         | Minimum UMI_3 count used in augmented MAF filter criteria                                                   | 1       |
| phase_id_col                | SAGE phase/ID column name                                                                                   | "LPS"   |
| phased_min_t_alt            | Minimum tumor alt count for phased variants                                                                 | 3       |


# Note about patient reports

Due to the uniqueness of the required analyses of different projects the `patient_reports.smk`  is not a one-click solution
for all projects. Therefore the workflow and code here will likely need to be adapted for your usage. Also included here is
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
As this code was developed on lymphoma patient samples, it can integrate the LymphGen classifier `lcr-modules` outputs. You can simply
set that to `false` in the config if you do not need it.

# Panel of Normals

To help exclude technical artifacts unique to a hybrid capture panel, it is recommended you sequence some cfDNA
from "healthy" individuals, and use their data to create a panel of normals. That VCF is required by the pipeline for filtering,
but if you don't have one you can provide an empty dummy file.

To assist in the creation of a PON, a workflow using GATK is included in this module: `make_somaticPON.smk`

# DLBCL Panel Hotspot and Notspot Lists

This pipeline was designed during a project using [DLBCL-STORM](https://github.com/morinlab/DLBCL-STORM "DLBCL-STORM GitHub repo"),
which has its own custom hybrid capture panel and hotspot/notspot lists.
