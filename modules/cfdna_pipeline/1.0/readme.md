# cfDNA Pipeline documentation

This pipeline coordinates the analysis of DNA sequencing data from hybrid capture assays of cfDNA. 

The pipeline consists of two main workflows, one aligning fastqs and deduplicating reads using UMIs,
and then variant calling and filtering using SAGE [https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md] as a caller.

The pipeline requires some resources not contained in this repo, links for some are provided below.

For GSC users these are probably floating around already.

The pipeline uses snakemake (v. 8.10.8) and conda.

# Resources 

## SAGE
https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md

If you need to download resources for SAGE check their downloads page:
https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_RESOURCES.md#variants


# Samplesheet requirements
The samplesheet has the following required columns

sample_id: must be unique
patient_id: must be the same as other samples from the same patient
tissue_status: tumor or normal
matched_normal: the sample_id of the matched normal sample. OR a value of "unmatched" if using the unmatched sample workflow.
genome_build: hg37 or hg38 (make sure you supply the right versions of resources in config)
seq_type: capture or genome (for lymphgen)
capture_space: name of panel, corresponds to dict in config

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
include: .../lcr-modules/modules/cfdna_pipeline/1.0/cfdna_UMI_workflow.smk
inlcude: .../lcr-modules/modules/cfdna_pipeline/1.0/cfDNA_SAGE_workflow.smk

# specify your rule all

rule all:
    input:
        rules.all_bams.input,
        rules.all_sage.input,
    threads: 1
    resources:
        mem_mb = 500

```
# Dependencies
The required dependencies are included as conda recipes in the envs directory, and snakemake
should automatically set them up for you. This should include any required python modules used
in any of the utils scripts.

# Note about Fastq inputs
You will see in the config that the pipeline requires an example path to your fastqs using a wildcard
in place of the sample id. This means you need to have all your fastqs nicely formatted in one input directory,
so if you need you can create an input directory for them and symlink them all there before running the pipeline.

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

.../lcr-modules/modules/cfdna_pipeline/1.0/patient_reports/report_template.ipynb
.../lcr-modules/modules/cfdna_pipeline/1.0/patient_reports.smk
.../lcr-modules/modules/cfdna_pipeline/1.0/patient_reports/compile_report.py

# Panel of Normals

The sage pipeline has a config item of "notlist", which can be any .vcf of genomic sites you want
to exclude from variant calling. It is suggested to creat this (or add to it) but creating a somatic
panel of normals using the included pipeline here, by feeding it matched normals.

The output vcf from that pipeline can be supplied as the "notlist", which result in common technical artifacts
and germline mutations being filterd from your data.
