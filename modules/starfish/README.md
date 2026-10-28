# starfish

# Purpose

The `starfish` is a Level 3 module that operates on `VCF` files to use Starfish and perform Aggregation for `capture, genome, mrna` data. It generates `VCF` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    starfish:
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

configfile: "../modules/starfish/2.0/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/starfish/2.0/starfish.smk"

rule all:
    input:
        rules._starfish_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/starfish/CHANGELOG.md)
