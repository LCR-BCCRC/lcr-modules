#!/bin/bash

# Launches a snakefile of your choice in test mode (for debugging)
# Usage: ./test.sh <snakefile.smk> <target_rule> "<snakemake_flags>"
# Example: ./test.sh example.smk example_all
# snakefile.smk The snakefile you want to run
# target_rule: The name of one of the target rules specified in one of the included Snakefiles
# snakemake_flags: One or more flags for the snakemake to run, specified inside quotation marks


# Default to all targets
snakefile=$1
TARGETS=${2:-all}
snakemake_flags=$3

./run.sh $snakefile $TARGETS  $snakemake_flags 

snakemake -np --dag -s $snakemake_flags $snakefile  $TARGETS| dot -Tpdf > dag_group.pdf 