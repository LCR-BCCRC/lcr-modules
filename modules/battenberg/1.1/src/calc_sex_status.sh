#!/usr/bin/env bash

# Infer sex of a patient from a normal genome bam using the ratio of X and Y chromosome reads
# If provided, the names of chrX and chrY will be used, otherwise they will be inferred from the header (this is not guaranteed to work)

set -euf -o pipefail

BAM="$1"
SAMPLE="${2:-UNKNOWN}"
X_CHROM="${3:-MISSING}"
Y_CHROM="${4:-MISSING}"
DEBUG="${5:-MISSING}"

if [[ $X_CHROM == "MISSING" ]]
then
    X_CHROM=$(samtools view -H ${BAM} |\
     sed -r 's/\S+:\S+/\n&/g' | perl -ne 's/\s+//g;print "$_\n"' | awk 'BEGIN{FS=":"} $1=="SN" && $2 ~ /X$/ {print $2}')
fi
if [[ $Y_CHROM == "MISSING" ]]
then
    Y_CHROM=$(samtools view -H ${BAM} |\
     sed -r 's/\S+:\S+/\n&/g' | perl -ne 's/\s+//g;print "$_\n"' | awk 'BEGIN{FS=":"} $1=="SN" && $2 ~ /Y$/ {print $2}')
fi
if [[ ! $DEBUG == "MISSING" ]]
then
    echo "DEBUG: x chromosome is named >$X_CHROM< and y chromosome is named >$Y_CHROM<"
fi


samtools idxstats $(readlink -e ${BAM}) \
    | awk -v sample="$SAMPLE" ' \
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
            ratio = chrY/chrX ; \
            if( ratio > 0.1) sex = "male"; \
            else sex = "female"; \
            print "sample", "chrX_count", "chrY_count", "sex"; \
            print sample, chrX, chrY, sex; \
        } \
    '
