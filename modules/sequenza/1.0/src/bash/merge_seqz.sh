#!/bin/bash

# Usage: 
#   merge_seqz.sh chr1.seqz.gz chr2.seqz.gz [...] | gzip > all.seqz.gz

set -euf

gzip -dc "$1" | head -1

zcat "$@" | egrep -v "^chromosome" 
