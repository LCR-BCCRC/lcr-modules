#!/bin/bash

# This script will filter out the non-standard chromosomes and non-ACTG characters from
# vcf files produced by lofreq. It accepts .qz-compressed files.
# Usage: 
#   lofreq_filter.sh input.vcf.gz | gzip > output.vcf.gz

INPUT_FILE="$1"

zcat "${INPUT_FILE}" \
  | awk '($4=="A" || $4 == "C" || $4=="T" || $4=="G" || /\#/)' \
  | perl -ne 'print if /^#|^(chr)*[\dX]+\s.+/' \
  | perl -ne 's/AF=/VAF=/g;s/ID=AF/ID=VAF/;print;' \
  | perl -F'\t' -lane '{
      $HEADER_DP = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">";
      $HEADER_DP4 = "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\">";
      $TO_HEADER = join "\t", "FORMAT", "TUMOR";
      $FORMAT_HEADER = "DP:DP4";
      ($DEPTH, $VAF, $SB, $DP4, $SOMATIC, $UQ) = split /;/, $F[7]; 
      ($DP, $DEPTH) = split /=/, $DEPTH;
      ($DP4_INFO, $DP4_FORMAT) = split /=/, $DP4;
      $FORMAT = join ":", $DEPTH, $DP4_FORMAT;
      if ($F[0] =~ /^#{2}/){
        print join "\t", @F; print join "\n", $HEADER_DP, $HEADER_DP4 if /HRUN/
      }elsif($F[0] =~ /^#{1}/){
        print join "\t", @F, $TO_HEADER
      }else{
        print join "\t", @F, $FORMAT_HEADER, $FORMAT
      }
      }'
