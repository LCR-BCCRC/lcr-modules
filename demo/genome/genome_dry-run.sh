#!/bin/bash

# Default to all targets
TARGETS=${@:-all}

snakemake --dryrun --cores 24 -s genome_Snakefile --printshellcmds --reason --use-conda $TARGETS
