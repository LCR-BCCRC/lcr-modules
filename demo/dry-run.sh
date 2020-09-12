#!/bin/bash

# Default to all targets
TARGETS=${@:-all}

snakemake --dryrun --cores 24 --printshellcmds --reason --use-conda $TARGETS
