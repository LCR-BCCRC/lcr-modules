#!/bin/bash

# This script will filter out the non-standard chromosomes and non-ACTG characters from
# vcf files produced by lofreq. It accepts .qz-compressed files, and outputs .gz files as well.
# Usage: 
#   lofreq_filter.sh input.vcf.gz | gzip > output.vcf.gz

INPUT_FILE="$1"

zcat "${INPUT_FILE}" \
	| awk '($4=="A" || $4 == "C" || $4=="T" || $4=="G" || /\#/)' \
	| grep -E -v "v1|v2|GL0"
