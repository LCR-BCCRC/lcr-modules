#!/usr/bin/env bash

set -o pipefail
set -e  # Exit immediately if a command fails

R --vanilla -q -e 'devtools::install_github(
    "Simon-Coetzee/motifBreakR",
    dependencies = FALSE,
    upgrade = "never"
)'

R --vanilla -q -e 'devtools::install_github(
    "morinlab/GAMBLR.data",
    dependencies = FALSE,
    upgrade = "never"
)'