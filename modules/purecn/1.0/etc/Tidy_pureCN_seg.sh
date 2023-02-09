#!/bin/bash

## instruction message
if [ "$1" == "-h" ] ; then
    echo -e "Using seg files generated from pureCN and converting them to a format for downstream reformatting.
Usage: `basename $0` -i <seg file> -o <output seg file>
    <seg file>: input seg file - type the path; 
    <out file>: output seg file - type the path"
    exit 0
fi

## command line arguments
seg=''
out=''
while [ $# -gt 0 ]
do
    case "$1" in
        -i) seg="$2"; shift;;
        -o) out="$2"; shift;;
    esac
    shift
done

echo -e "Running... $seg"
echo -e "ID\tchrom\tstart\tend\tnum.mark\tseg.mean" > $out

echo $seg
awk '(NR>1)' $seg | awk '{print $1"\tchr"$2"\t"int($3)"\t"int($4)"\t"$5"\t"log($7)/log(2)-1}' | sed 's/chr23/chrX/g' | sed 's/chr24/chrY/g' | sed 's/chrchr/chr/g' | sed 's/-inf/-2/g' | sort -k1,1 -k2,2 -k3,3n -k4,4n -V | awk '($4 > $3)' |  awk '$3+1!=$4 {print}' >> $out

echo -e "Finished writing... $out"

