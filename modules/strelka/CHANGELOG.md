# Changelog

All notable changes to the `strelka` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0] - 2020-06-03

This release was authored by Helena Winata.

- `_strelka_run` outputs to a directory because strelka creates vcf files in `${STRELKA_ANALYSIS_PATH}/results/variants`
- `_strelka_dispatch` is used to call outputs based on the `{pair_status}` wildcard since outputs files are automatically generate.
    - somatic workflow outputs `somatic.snvs.vcf.gz` and `somatic.indels.vcf.gz`
    - germline workflow outputs `variants.vcf.gz`
- It is recommended to run Strelka with [candidatesSmallIndels.vcf] (https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#capabilities) from Manta variant caller.
    - Specify the `candidateSmallIndels.vcf` file under `CFG["inputs"]["candidate_small_indel"]`
