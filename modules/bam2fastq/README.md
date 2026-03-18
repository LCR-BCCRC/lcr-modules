# Purpose

The `bam2fastq` is a Level 1 module that operates on `bam` or `cram` files to use picard tools and generate `fastq` files as outputs.

# Example

To run this module, have config and snakefile in the current directory. The example config:

```yaml
lcr-modules:
    _shared:
        lcr-modules: "../" # path to the lcr-modules directory relative to the location of the running snakefile
        lcr-scripts: "../../lcr-scripts/" # path to the lcr-scripts directory relative to the location of the running snakefile
        root_output_dir: "results/" # path to the location of the outputs
        scratch_directory: "scratch/" # path to the location of the intermediate files

    bam2fastq:
        inputs:
            sample_bam: "data/{sample_id}.bam" # path to the input files
        temp_outputs: True # fastq outputs will be temporary
```

The example snakefile:

```python
#!/usr/bin/env snakemake


##### SETUP #####
import oncopipe as op
# define sample table
SAMPLES = op.load_samples("data/samples.tsv")


##### REFERENCE_FILES WORKFLOW #####
subworkflow reference_files:
    workdir:
        "reference/" # path to the directory where all the references will be stored
    snakefile:
        "../workflows/reference_files/2.4/reference_files.smk"
    configfile:
        "../workflows/reference_files/2.4/config/default.yaml"


##### CONFIGURATION FILES #####
# Load module-specific configuration
configfile: "../modules/bam2fastq/1.2/config/default.yaml" # path to the default config relative to the location of the running snakefile
# Load project-specific config from the example above
configfile: "capture.yaml"


##### CONFIGURATION UPDATES #####
# Use all samples as a default sample list for each module
config["lcr-modules"]["_shared"]["samples"] = SAMPLES

##### MODULE SNAKEFILES #####
# Load module-specific snakefiles
include: "../modules/bam2fastq/1.2/bam2fastq.smk" # path to the module snakefile relative to the location of the running snakefile

##### TARGETS ######
rule all:
    input:
        rules._bam2fastq_all.input # request the output files
```
