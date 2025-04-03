#!env bash
set -o pipefail
set -e # Exit on error

gamblr_sha="d8d643face5b67bdfe850aabcbfe26f82254a825"

R --vanilla -q -e \
    'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR", ref="d8d643face5b67bdfe850aabcbfe26f82254a825", dependencies = TRUE, upgrade = "never"); devtools::install_github("morinlab/ggsci", dependencies = TRUE, upgrade = "never"); if(!require(GAMBLR)) stop("GAMBLR is not installed")'
