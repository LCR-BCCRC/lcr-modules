# Changelog

All notable changes to the `purecn` module will be documented in this file.

## [1.0] - 2026-06-08

This release was authored by translating the CCSRI_1500 exomes PureCN pipeline
(`/projects/dscott_prj/CCSRI_1500/exomes/snakefile_purecn.smk`, J. Wong) into a
standalone lcr-module.

- Estimates tumour purity, ploidy, integer copy number and LOH with PureCN
  (bioconductor-purecn 2.8.1), tumour-only.
- Consumes CNVkit output from the `cnvkit/1.0` module (`.cnr` + BAF `call.cns`).
- Consumes the panel of normals from the `panel_of_normals/1.0` module
  (intervals, mapping_bias.rds, denovo normalDB.rds, Mutect2 PON VCF). The PON
  construction was split out of this module so cnvkit and purecn share one
  normal set.
- Runs Mutect2 on the tumour (vs the PON) for BAF, then PureCN in two modes per
  tumour (cnvkit mode + denovo mode).
- Parametrized for grch37 and hg38 (genome_builds_map -> hg19 / hg38); all
  references resolved through the reference_files workflow.
- Reuses the `mutect2/2.0` gatk/bcftools envs and the `cnvkit/1.0` env.
