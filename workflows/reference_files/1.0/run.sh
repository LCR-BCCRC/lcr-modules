#!/bin/bash

snakemake --dryrun --debug --cores 1 --printshellcmds --reason --use-conda --snakefile prepare_reference_files.smk --config reference_directory=/projects/bgrande/reference_files
