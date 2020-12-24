#!/usr/bin/env bash

# Infer sex of a patient from a normal genome bam using the ratio of X and Y chromosome reads
# If provided, the names of chrX and chrY will be used, otherwise they will be inferred from the header (this is not guaranteed to work)

set -euf -o pipefail

BAM="$1"
REF="$2"
SAMPLE="${3:-UNKNOWN}"
X_CHROM="${4:-MISSING}"
Y_CHROM="${5:-MISSING}"
DEBUG="${6:-MISSING}"

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


X_READS=$(samtools view -@ 8 -T $REF $BAM $X_CHROM | wc -l)
Y_READS=$(samtools view -@ 8 -T $REF $BAM $Y_CHROM | wc -l)


ratio=$((100 * $Y_READS/$X_READS))
sex="female"
if [[ $ratio -gt 10 ]]
then
    sex="male"
fi
printf "sample\tchrX_count\tchrY_count\tsex\n"
printf "$SAMPLE\t$X_READS\t$Y_READS\t$sex\n"
