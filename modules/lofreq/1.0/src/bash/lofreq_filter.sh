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
  | perl -ne '
              # Add 2 new rows to the description and 2 new columns in the header
              if(/^#/){
                if(/##INFO=<ID=DP,.+\n/){
                  $DP=$&;
                };
                $DP =~ s/INFO/FORMAT/;
                print $DP if /min_dp/;
                if(/##INFO=<ID=DP4,.+\n/){
                  $DP4=$&;
                };
                $DP4 =~ s/INFO/FORMAT/;
                print $DP4 if /min_dp/;
                if(/^#CHROM.+/){
                  s/$&/$&\tFORMAT\tTUMOR/;
                };
                print;
              }
              # For each feature, add FORMAT column with descriptors and populate TUMOUR column with depth, reads counts
              else{
                my @data = map { chomp; [ split /=|;/ ] } $_;
                $NEW_ROW = "$_\tDP:DP4\t$data[0][1]:$data[0][7]\n";
                print $NEW_ROW;
              }'
