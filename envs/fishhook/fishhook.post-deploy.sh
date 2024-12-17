#!env bash
set -o pipefail
set -e # Exit on error

R --vanilla -q -e 'devtools::install_github("mskilab/gUtils", dependencies = FALSE, upgrade = "never")'
R --vanilla -q -e 'devtools::install_github("mskilab/gTrack", dependencies = FALSE, upgrade = "never")'
R --vanilla -q -e 'devtools::install_github("mskilab/fishHook", dependencies = FALSE, upgrade = "never")'