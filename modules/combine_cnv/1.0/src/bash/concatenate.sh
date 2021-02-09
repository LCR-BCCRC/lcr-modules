#!/usr/bin/env bash

# This script will concatenate all filtered samples in specified directory and output them in a single .seg file.
# Usage: 
# bash concatenate.sh <input_directory> <output_file>

INPUT_DIR="$1"
SEG_PATH="$2"
counter=0

echo "Finding filtered seg files in..."

for i in $(find ${INPUT_DIR} -type f -name "*.filtered.seg");
do
counter=$(($counter +1))
if [ $counter = 1 ]
then
cat $i | head -1 > $SEG_PATH
fi

cat $i | grep -v "sample" >> $SEG_PATH


done

echo "Successfully merged $counter sample(s)."
