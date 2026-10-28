# picard_qc

# Purpose

The `picard_qc` is a Level 1 module that operates on `BAM/CRAM` files to use Picard tools and perform QC for `capture, genome, mrna` data. It generates `TSV` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    picard_qc:
        inputs:
            sample_bam: "data/{sample_id}.bam"
```

The example snakefile:

```python
#!/usr/bin/env snakemake

import oncopipe as op

SAMPLES = op.load_samples("data/samples.tsv")

subworkflow reference_files:
    workdir:
        "reference/"
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"

configfile: "../modules/picard_qc/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/picard_qc/1.0/picard_qc.smk"

rule all:
    input:
        rules._picard_qc_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/picard_qc/CHANGELOG.md)
