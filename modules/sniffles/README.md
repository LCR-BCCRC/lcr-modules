# sniffles

# Purpose

The `sniffles` is a Level 2 module that operates on `BAM/CRAM` files to use Sniffles and perform Structural variants long reads for `promethION` data. It generates `VCF` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    sniffles:
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

configfile: "../modules/sniffles/1.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/sniffles/1.0/sniffles.smk"

rule all:
    input:
        rules._sniffles_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/sniffles/CHANGELOG.md)
