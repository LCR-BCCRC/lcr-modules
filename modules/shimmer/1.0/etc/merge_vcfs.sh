#!/bin/bash

set -euf -o pipefail
LC_ALL=C

# Output header from one file
zgrep "^#" $1

# Contatenate the non-header portions of all files
zgrep -v -h "^#" $@
