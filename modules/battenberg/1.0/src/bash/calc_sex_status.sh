#!/usr/bin/env bash

set -euf -o pipefail

BAM="$1"
X_CHROM="${2:-X}"
Y_CHROM="${3:-Y}"


samtools idxstats $(readlink -e ${BAM}) \
    | awk -v sample="SAMPLE" ' \
        BEGIN { \
            FS=OFS="\t" \
        } \
        $1 == "'${X_CHROM}'" { \
            chrX = $3 \
        } \
        $1 == "'${Y_CHROM}'" { \
            chrY = $3 \
        } \
        END { \
            print "sample", "chrX_count", "chrY_count"; \
            print sample, chrX, chrY; \
        } \
    '
