#!/bin/bash

set -euf -o pipefail
LC_ALL=C

# Output header from one file
zgrep "^Index" $1

# Contatenate the non-header portions of all files
zgrep -v -h "^Index" $@
