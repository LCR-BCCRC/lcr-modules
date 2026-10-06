# Changelog

All notable changes to the `ega_download` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0] - 2026-10-05

This release was authored by Laura Hilton.

### Breaking Changes

- **pyega3 replaced by SFTP + crypt4gh.** The new EGA Live Outbox distribution method
  delivers crypt4gh-encrypted files via SFTP. This module no longer uses pyega3.
- **No `credentials_file`.** Authentication for SFTP download uses your registered SSH
  public key (`id_rsa.pub` or `id_ed25519.pub`), which must be shared with EGA to gain
  Live Outbox access. `private_key` is used exclusively for crypt4gh decryption.
- **Sample metadata loaded via `op.load_samples()`**, not a module-managed CSV. The
  sample table must include a `file_name` column containing the EGA file name as it
  appears on the SFTP outbox (e.g. `LBCL_scAtlas_2117-01_1_WES-tu_R1.fastq.gz`),
  without the `.c4gh` suffix. The required columns are: `sample_id`, `patient_id`,
  `tissue_status`, `seq_type`, `file_name`. One row per physical file.
- **`egad` is now a config value**, not a sample table column. Set it in your project
  config (e.g. `egad: EGAD50000001491`).
- **Output files are named by `{seq_type}/{file_name}`**, preserving the EGA file name.
  v1.0 outputs used `{seq_type}/{study_id}/{sample_id}{read}.{file_format}`.

### New Features

- **True parallelism.** Each file is an independent Snakemake job. With `snakemake
  --jobs N`, up to N files download simultaneously. v1.0 pyega3 ran serially.
- **Two-step resilient pipeline.** Download (`.c4gh`) and decrypt are separate rules.
  If decryption fails, re-running only re-decrypts without re-downloading.
- **Supports new EGA accession format** (EGAN50..., EGAF50...) from EGA Live Outbox.
- **New seq_types in pairing_config**: `snrna`, `snatac`, `lpwgs` added alongside the
  existing `genome`, `capture`, and `mrna`.

### Migration from v1.0

1. Obtain EGA Live Outbox SFTP access for your dataset.
   Register your SSH public key with EGA (`ssh-keygen -t ed25519` if needed).
2. Generate a crypt4gh key pair and submit the public key to EGA for re-encryption:
   `crypt4gh-keygen --sk private.key --pk public.key`
3. Build a sample metadata TSV with columns: `sample_id`, `patient_id`, `tissue_status`,
   `seq_type`, `file_name`. The `file_name` value for each row is the file name as listed
   in the EGA outbox or EGA portal metadata (without the `.c4gh` suffix).
4. In your Snakefile, load samples with `op.load_samples()` and assign to
   `config["lcr-modules"]["_shared"]["samples"]` as for any other module.
5. In your config, set `egad`, `private_key`, and `sftp_server` instead of the v1.0
   `credentials_file`, `study_id`, and `out_file_type` keys.

### Example Sample Metadata TSV (v2.0)

| sample_id         | patient_id | tissue_status | seq_type | file_name                                   |
|-------------------|------------|---------------|----------|---------------------------------------------|
| 2117-01_1_WES-tu  | 2117-01    | tumour        | capture  | LBCL_scAtlas_2117-01_1_WES-tu_R1.fastq.gz  |
| 2117-01_1_WES-tu  | 2117-01    | tumour        | capture  | LBCL_scAtlas_2117-01_1_WES-tu_R2.fastq.gz  |

---

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
