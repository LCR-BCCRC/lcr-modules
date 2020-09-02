#!/bin/bash

#set -euf -o pipefail
#LC_ALL=C

# Output header from one file
grep "start" $1

# Contatenate the non-header portions of all files
grep -v -h "start" $@
