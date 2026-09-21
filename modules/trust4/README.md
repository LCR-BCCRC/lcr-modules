# trust4

# Purpose

The `trust4` module is a Level 2 module that operates on `BAM` files to use TRUST4 ([liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4)) and reconstruct BCR/TCR immune repertoires (V/D/J/C genes, CDR3 sequences) from bulk RNA-seq for `mrna` data. It generates `TSV`/`FASTA` files as outputs.

This is a **bulk RNA-seq v1**: it takes an existing coordinate-sorted, indexed BAM as input (by default, `modules/star/1.4`'s own final output) and does not run TRUST4's own barcode/UMI/10x single-cell machinery, nor its raw-FASTQ input mode. TRUST4 operates on one BAM per `sample_id` -- there is no tumour/normal pairing concept anywhere in this module.

TRUST4 is bioconda-installable (`trust4=1.1.10`) with no licensing restriction, unlike `modules/mhc_hammer/1.0`'s Novoalign/HLA-HD/VEP dependencies -- this module installs it via a normal conda/container environment, no user-supplied external script directory needed.

## Reference files

TRUST4's own GitHub repo ships small, pre-built, plain-text reference FASTAs directly in-repo (not via bioconda): a genome-coordinate-aware V/D/J/C gene FASTA per genome build (`hg19_bcrtcr.fa`/`hg38_bcrtcr.fa`), and a species-only IMGT reference (`human_IMGT+C.fa`). This module downloads both automatically, pinned to a specific commit (`options.trust4_repo_commit`) for reproducibility -- no manual reference-building step needed for standard human genome builds (`grch37`/`hg19`/`grch38`/`hg38`, mapped via `options.trust4_bcrtcr_map`).

## What's not included in v1

- Raw FASTQ input (TRUST4's own `-1/-2/-u` alternative to `-b`) -- this module is BAM-only.
- Barcode/UMI/10x single-cell support (`--barcode`, `--barcodeWhitelist`, `--UMI`, and the resulting `trust-barcoderep.pl`-produced `_barcode_report.tsv`/`_barcode_airr.tsv` outputs) -- bulk RNA-seq only.
- Mouse genome builds (TRUST4 also ships `GRCm38_bcrtcr.fa`/`GRCm39_bcrtcr.fa`/`mouse_IMGT+C.fa`) -- a trivial addition to `options.trust4_bcrtcr_map` if ever needed, not wired up now.
- `--assembleWithRef`/`--stage`/`--outputReadAssignment` and other advanced TRUST4 options -- not exposed as module options; open an issue/PR if you need one.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    trust4:
        inputs:
            # Defaults to modules/star/1.4's own output -- override only if your own
            # deployment uses a different star version/path, or a different RNA-seq aligner
            sample_bam: "star-1.4/99-outputs/bam/{seq_type}--{genome_build}/{sample_id}.bam"
            sample_bai: "star-1.4/99-outputs/bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

configfile: "../modules/trust4/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/trust4/1.0/trust4.smk"

rule all:
    input:
        rules._trust4_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/trust4/CHANGELOG.md)
