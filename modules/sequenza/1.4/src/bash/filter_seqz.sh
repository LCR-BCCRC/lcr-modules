#!/bin/bash

# Usage: 
#   merge_seqz.sh all.seqz.gz dbsnp.pos [rare_vars.pos] | gzip > all.filt.seqz.gz

SEQZ_FILE="$1"
DBSNP_POS_FILE="$2"
RARE_VARIANTS_TMP=$(mktemp /tmp/merge_seqz.sh.XXXXXX)
RARE_VARIANTS="${3:-${RARE_VARIANTS_TMP}}"

BUFFER_SIZE="${BUFFER_SIZE:-20G}"
export LC_ALL=C  # Set sort locale to be consistent
# Check the DPSNP_POS_FILE is sorted
sort -c $2 || { echo "Provided DBSNP_POS file is not sorted. Check that your \'LC_ALL\' envionmental variable is set to \'$LC_ALL\'" && exit 1; }

set -euf -o pipefail

zcat "${SEQZ_FILE}" \
	| egrep -v "^chromosome" \
	| awk 'BEGIN {FS="\t"} $9 == "het" {print $1 ":" $2}' \
	| sort -S "${BUFFER_SIZE}" --parallel 10 \
	| comm - "${DBSNP_POS_FILE}" -2 -3 \
	| tr ":" "\t" \
	> "${RARE_VARIANTS}"

zgrep -v -F -f "${RARE_VARIANTS}" "${SEQZ_FILE}"

rm -f "${RARE_VARIANTS_TMP}"
