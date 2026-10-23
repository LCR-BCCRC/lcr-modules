#!env bash
set -o pipefail
set -e # Exit on error

GAMBLR_HELPERS_REF="${GAMBLR_HELPERS_REF:-fdde14c}"

R --vanilla -q -e "remotes::install_github(
    'morinlab/GAMBLR.helpers',
    ref = '${GAMBLR_HELPERS_REF}',
    dependencies = NA,
    upgrade = 'never'
)"

