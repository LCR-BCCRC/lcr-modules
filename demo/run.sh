#!/bin/bash

# Number of cores set by first command-line argument (default is 4)
NUM_CORES=${1:-4}

snakemake --cores "${NUM_CORES}" --printshellcmds --use-conda all
