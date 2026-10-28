# ecotyper

# Purpose

The `ecotyper` is a Level 3 module that operates on `TSV` files to use EcoTyper (Steen et al, 2021) and perform Microenvironment analysis for `mrna` data. It generates `TSV` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    ecotyper:
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

configfile: "../modules/ecotyper/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/ecotyper/1.0/ecotyper.smk"

rule all:
    input:
        rules._ecotyper_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/ecotyper/CHANGELOG.md)
