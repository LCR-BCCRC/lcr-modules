#!/usr/bin/env bash
set -o pipefail
set -e # Exit on error

R --vanilla -q -e 'remotes::install_github("caravagnalab/CNAqc", dependencies = TRUE, upgrade = "never")'
