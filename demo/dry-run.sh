#!/bin/bash

snakemake --dryrun --cores 1 --printshellcmds --reason --use-conda all
