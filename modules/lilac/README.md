# lilac

# Purpose

The `lilac` module is a Level 2 module that operates on `BAM` files to use HMF Tools' LILAC ([hartwigmedical/hmftools](https://github.com/hartwigmedical/hmftools/tree/master/lilac)) for HLA class I typing (`A`/`B`/`C`, 4-digit resolution), somatic mutation assignment per allele, allelic imbalance, and complete allele loss detection, for `capture, genome` (WES/WGS) data. It generates `TSV`/`VCF` files as outputs.

This is the WGS-native alternative to `modules/mhc_hammer/1.0`: LILAC is explicitly tested and validated at 30-40x germline / 100x tumour WGS depth, whereas MHC Hammer's own design was validated and tuned for WES and has repeatedly needed depth-filter tuning on WGS data in practice (see `modules/mhc_hammer/1.0/README.md`).

Unlike `mhc_hammer`, LILAC does its own HLA-region fragment filtering directly against a full (unsliced) BAM -- there is no personalised-reference construction, realignment, or allele-specific BAM splitting step. A single `lilac` invocation per matched tumour/normal pair does HLA typing, and (when the optional inputs below are provided) tumour copy number and somatic mutation assignment, all at once.

LILAC is distributed under the **GPLv3 licence** -- unlike `mhc_hammer`'s upstream (Cancer Research Horizons academic-use licence), there is no redistribution/modification restriction, no user-supplied-clone requirement. It's installed as a normal bioconda package (`hmftools-lilac`), like any other tool this repo wraps.

## Prerequisites (read before use)

Nothing needs to be supplied manually. `_lilac_download_reference` fetches LILAC's own 3 required resource files (`lilac_allele_frequencies.csv`, `hla_ref_nucleotide_sequences.csv`, `hla_ref_aminoacid_sequences.csv`, ~30MB combined) from `www.bcgsc.ca/downloads/morinlab/hmftools-references/lilac/` -- the same pre-extracted, per-tool reference mirror `modules/hmftools/1.1` and `modules/sage/1.1` already use for this upstream tool suite. One cohort-wide download, not per genome build (all 3 files are build-independent). The reference genome FASTA and the `lilac` binary itself are handled the standard lcr-modules way (`reference_files()`, `conda_envs`).

Every patient in your sample table must have **exactly one** germline WES/WGS sample (`tissue_status: normal`) per `patient_id`/`seq_type`/`genome_build` combination that you want typed -- oncopipe's own pairing (`CFG["paired_runs"]`, narrowed to `pair_status == "matched"` only) handles this; a tumour sample without a real matched germline in the same patient is never processed, since HLA typing needs the patient's own germline sample.

## Optional inputs

Both are optional per pair -- missing either one just means LILAC runs without the corresponding output (HLA typing itself is unaffected either way):

- `inputs.gene_copy_number`: a gene-level copy number TSV (minimum copy number and minimum minor allele copy number per gene), e.g. PURPLE's own `cnv.gene.tsv` output (`modules/hmftools/1.1`). Enables tumour allele-specific copy number output.
- `inputs.somatic_vcf`: a somatic SNV/indel VCF with `TUMOR`/`NORMAL`-labelled samples, e.g. SAGE's own combined VCF (`modules/sage/1.1`) or PURPLE's. Enables per-allele somatic mutation assignment (`SomaticMissense`, `SomaticNonsenseOrFrameshift`, etc. in the output TSV, plus an annotated `<tumour_id>.lilac.somatic.vcf.gz`).

## What's not included in v1

- LILAC's own native tumour-only and germline-only run modes -- this module always runs paired mode (germline + tumour BAM together), matching `mhc_hammer`'s own v1 simplification, even though LILAC itself could type a germline-only or tumour-only sample. A real, low-cost future extension, not implemented yet.
- BAM slicing/subsetting to the HLA region before running LILAC -- not needed for correctness (LILAC does its own internal region filtering against a full BAM), only ever worth adding as a runtime optimisation if it becomes a real bottleneck in practice.
- A cohort-wide combined table across patients, and a VCF-to-MAF conversion of the somatic VCF output -- both should be much simpler here than the equivalent `mhc_hammer` features (LILAC's own per-sample TSV is already close to cohort-table-shaped, and its somatic VCF already has real genomic coordinates since it's just annotating an existing real-coordinate VCF, not calling against a personalised allele-specific reference) -- worth adding as a fast follow-up, not required for v1.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    lilac:
        inputs:
            sample_bam: "data/{sample_id}.bam"
            sample_bai: "data/{sample_id}.bam.bai"
            # Optional -- omit or leave as "" to run without tumour copy number / somatic mutation output
            gene_copy_number: "results/hmftools-1.1/99-outputs/purple_output/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purple.cnv.gene.tsv"
            somatic_vcf: "results/sage-1.1/99-outputs/combined/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sage.combined.vcf.gz"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/lilac/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/lilac/1.0/lilac.smk"

rule all:
    input:
        rules._lilac_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/lilac/CHANGELOG.md)
