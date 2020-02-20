# Morin Lab Snakemake Pipelines

## Setup

These setup instructions assume you already have `snakemake` installed and have a Snakefile for your project.

1) Install the `pipeline_utils`modules. 

```bash
pip install -e path/to/pipelines/pipeline_utils
```

2) Add the path of the pipelines repository/directory to your configuration. 

```yaml
pipelines:
    repository: "path/to/pipelines/"
```

3) Load the pipeline configuration first and then include the pipeline Snakefile. 

```yaml
configfile: "src/pipelines/strelka/1.0/config.yaml"
include:    "src/pipelines/strelka/1.0/Snakefile"
```