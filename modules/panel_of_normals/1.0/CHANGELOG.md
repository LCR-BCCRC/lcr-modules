# Changelog

## [1.0] - 2026-06-08

PureCN-side panel of normals, split out of the purecn/1.0 module so cnvkit and
purecn draw the PON from one normal set. Built per
{seq_type}--{genome_build}/{capture_space} from every normal in the samples table
(e.g. normal exomes injected from gambl_metadata_normals.tsv). Emits baits
intervals, mapping_bias.rds, denovo normalDB.rds, and a Mutect2 PON VCF.
Translated from /projects/dscott_prj/CCSRI_1500/exomes (J. Wong).
