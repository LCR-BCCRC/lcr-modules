# TEMPEST

TEMPEST (Temporal Evaluation of Mutations in Plasma cfDNA with Error-suppressed Somatic Tracking) is a pipeline designed
to process hybrid capture sequencing samples from liquid biopsies of cancer patients.

## Changelog:

### v 2.0.0

- reworked UMI deduplication using fgbio's suggested workflow, as a result the pipeline is generally more memory efficient and
can handle larger input fastqs, like from WES. 
- removed custom base pair flagging scripts that were slow in favour of fgbio's more optimized solutions.
- the fastq and variant calling workflows are now more modular, and any input bams can fed into the variant calling workflow.
- bumped the SAGE version to v.4.1

### v 1.6.0

- introduced the ability to filter variants based on background mutation rates, if an index is provided. 
- made different artifact filtering methods optional.