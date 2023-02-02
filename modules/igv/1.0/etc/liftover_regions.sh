#!/bin/bash

# Use CrossMap.py to convert genomic coordinates between GRCh37 and GRCh38 for MAF, VCF, BED or BEDPE files.
# 
# Usage: liftover_regions.sh \
# <input_regions> \
# <input_type> \
# <run_build> \
# <output_file> \
# <chain_file> \
# <target_ref> \
# <target_build>

input_regions=$1
input_type=$2
regions_build=$3
target_build=$4
output_file=$5
chain_file=$6
target_ref=$7

echo "Input regions file:       $input_regions"
echo "Input regions type:       $input_type"
echo "Input regions build:      $regions_build"
echo "Target genome build:      $target_build"
echo "Output file:              $output_file"
echo "Chain file:               $chain_file"
echo "Target reference:         $target_ref"

intermediate_output_file=$(echo $output_file)_int

# MAFs
# Check genome build of incoming MAF file to determine what build it needs to be changed to
if [ "$input_type" == "maf" ] ;
then
    echo "Proceeding with MAF input..."
    if [ $regions_build == $target_build ] ;
    then
        echo "WARNING: Input regions file $input_regions is already $target_build. Copying contents of $input_regions to $output_file";
        cut -f 1,5,6,7,9,10,11,13,16 $input_regions > $output_file 
    else
        echo "Input regions file $input_regions does not appear to be $target_build. Proceeding with conversion to $target_build"
        echo "CrossMap.py maf $chain_file $input_regions $target_ref $target_build $output_file"
        CrossMap.py maf $chain_file $input_regions $target_ref $target_build $intermediate_output_file
        cut -f 1,5,6,7,9,10,11,13,16 $intermediate_output_file > $output_file
        rm $intermediate_output_file
    fi
    echo "Finished MAF block."
fi

if [ "$input_type" == "bed" ] ;
    then
        echo "Proceeding with BED input..."
        if [ $regions_build == $target_build ] ;
        then
            echo "WARNING: Input regions file $input_regions is already $target_build. Copying contents of $input_regions to $output_file";
            cat $input_regions > $output_file
        else
            echo "Input regions file $input_regions does not appear to be $target_build. Proceeding with conversion to $target_build"
            echo "CrossMap.py bed $chain_file $input_regions $output_file"
            CrossMap.py bed $chain_file $input_regions $output_file
        fi
        echo "Finished BED block."
fi      

echo "End of bash script"
