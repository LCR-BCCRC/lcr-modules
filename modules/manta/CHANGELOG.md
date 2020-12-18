# Changelog

All notable changes to the `manta` module will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.31] - 2020-12-17

This minor version increase enables processing of cram files. The input is expected to be named ".bam" and ".bai" (or symlinks with that naming pointing to the cram and crai). 

## [2.2] - 2020-07-16

This release was authored by Bruno Grande.

- All `rules` references were wrapped with `str()`.

## [2.1] - 2020-06-06

This release was authored by Bruno Grande.

- The default configuration file was updated to include '__UPDATE__'.

## [2.0] - 2020-05-26

This release was authored by Bruno Grande.

- The module was updated to follow the latest best practices for `lcr-modules`. For instance, the "switches" (_i.e._ the dictionaries provided to the `switch_on_wildcard()` function) are stored in the configuration file. Another change is the use of the `reference_files` workflow.

- The augment_manta_vcf.py script is used instead of the older calc_manta_vaf.py script, which was written with the `somaticSV` VCF files in mind. The newer script is aware of which samples (_i.e._ tumour and/or normal) are present in each Manta output VCF file (_e.g._ `diploidSV`, `somaticSV`). This awareness also makes it simpler for this script to update the sample IDs (rather than relying on an awk command).

- The Snakemake checkpoint was dropped in version 2.0. Instead, the input file function for the `manta_dispatch` rule predicts which the Manta output VCF files based on the run. For example, if the `--rna` option is used, then the module will expect a `rnaSV` output file. For more details, check out the `_manta_predict_output` function in the module.

## [1.0] - 2020-04-22

This release was authored by Bruno Grande.

- A BED file consisting of the entire main chromosomes is provided to speed up computation as per [this recommendation](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#improving-runtime-for-references-with-many-short-contigs-such-as-grch38). As per the Manta user guide, the BED file is bgzip-compressed and tabix-indexed.

- The `manta` module supports both `paired` and `unpaired` analyses.

- To avoid duplicating the `_manta_configure` rule for paired and unpaired analyses, the module dynamically generates the shell command using the `switch_on_*()` functions. Namely, the `normal_bam_arg` and `tumour_bam_arg` params are determined based on the value of the `pair_status` wildcard and the `seq_type` wildcard, respectively. The values that will get loaded are specified in `default.yaml` under `switches`.

- In an earlier version of this module, the `normal_bam` input file would be omitted if `pair_status` was `no_normal` using `switch_on_wildcard()`. However, it was clunky because it had to use `unpack()` in case it returned no file. As a workaround, since the argument for the normal BAM file is omitted if the run is `no_normal` (see note on `normal_bam_arg` above), I added the `_manta_input_bam_none` rule to create an empty file to make sure Snakemake doesn't complain about an missing input file. I'm not happy with this approach, but it's cleaner than the messy code below.

  ```python
  unpack(op.switch_on_wildcard("pair_status", {
      "_default" : {"normal_bam": CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"},
      "no_normal" : {}
  }))
  ```

- For [RNA-seq](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#rna-seq) and [capture sequencing](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#exometargeted) data, Manta is run in [high-sensitivity mode](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#high-sensitivity-calling) as per the documentation.

- The `_manta_fix_vcf_ids` rule is there to ensure that the sample IDs in the VCF header match the sample IDs in Snakemake. This is necessary because Manta extracts the sample IDs from the read groups in the input BAM files, which may not match the sample IDs provided to Snakemake.

- The `_manta_calc_vaf` rule calculates the variant allele fraction (VAF) based on the reference and alternate allele counts in the VCF files. This rule uses the `calc_manta_vaf.py` script in `lcr-scripts`. The VAF information is appended to the `INFO` column, which is available in the BEDPE files. The output of this rule is what is symlinked into `99-outputs/vcf/`.

- The `manta` module uses a checkpoint. This was necessary because the output of Manta depends on how it's run. Specifically, the VCF files it produces depends on whether it was run in paired more and whether it was run with `--rna`. Accordingly, the checkpoint asks for the final VCF and BEDPE output files for the ones that are actually created. The function behind this, `_get_manta_files()`, also omits any empty VCF files from the BEDPE creation.
