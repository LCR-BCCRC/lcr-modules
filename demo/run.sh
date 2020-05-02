#!/bin/bash

snakemake --dryrun --cores 2 --printshellcmds --reason --use-conda _manta_all
