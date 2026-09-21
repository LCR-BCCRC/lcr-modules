# Changelog

All notable changes to the `ega_download` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2021-10-06

This release was authored by Kostiantyn Dreval.

- The input is a sample table of the following example:
| patient_id | sample_id       | file_name            | file_format | seq_type | genome_build | EGAS           | EGAD           | EGAN           | EGAF           |
|------------|------------------|----------------------|-------------|----------|--------------|----------------|----------------|----------------|----------------|
| E003       | E003.TI.FFPE.WES | E003_TI_FFPE_WES_R1  | fq.gz       | capture  | grch37       | EGAS00001006927 | EGAD00001011369 | EGAN00004220717 | EGAF00007894390 |
| E003       | E003.TI.FFPE.WES | E003_TI_FFPE_WES_R2  | fq.gz       | capture  | grch37       | EGAS00001006927 | EGAD00001011369 | EGAN00004220717 | EGAF00007894391 |
| E004       | E004.TI.FFPE.WES | E004_TI_FFPE_WES_R1  | fq.gz       | capture  | grch37       | EGAS00001006927 | EGAD00001011369 | EGAN00004220718 | EGAF00007894392 |
| E004       | E004.TI.FFPE.WES | E004_TI_FFPE_WES_R2  | fq.gz       | capture  | grch37       | EGAS00001006927 | EGAD00001011369 | EGAN00004220718 | EGAF00007894393 |
| E006       | E006.TI.FFPE.WES | E006_TI_FFPE_WES_R1  | fq.gz       | capture  | grch37       | EGAS00001006927 | EGAD00001011369 | EGAN00004220719 | EGAF00007894394 |


- The meaning of EGA id columns is to allow reproducible and unambiguous interpretation
of which sample is stored in the EGA under which accesion ID.
- The file_name value in the table above can be retrieved in the output of `pyega3 -cf </Path/To/CREDENTIALS_FILE> files EGAD<NUM>`. This will specify upfront the expected naming of
the sample after download from EGA.
