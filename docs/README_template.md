# {module_name}

# Purpose

The `{module_name}` is a Level {level} module that operates on `{input_format}` files to use {tool_name} and perform {module_purpose} for `{seq_types}` data. It generates `{output_format}` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../"
        lcr-scripts: "../../lcr-scripts/"
        root_output_dir: "results/"
        scratch_directory: "scratch/"

    {module_name}:
        inputs:
            sample_bam: "data/{{sample_id}}.bam"
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

configfile: "../modules/{module_name}/{version}/config/default.yaml"
configfile: "my_config.yaml" # the path to config file from the previous example

config["lcr-modules"]["_shared"]["samples"] = SAMPLES

include: "../modules/{module_name}/{version}/{module_name}.smk"

rule all:
    input:
        rules._{module_name}_all.input
```

# Changelog

See the full changelog [here](https://github.com/LCR-BCCRC/lcr-modules/blob/master/modules/{module_name}/CHANGELOG.md)
