#!env bash
set -o pipefail
set -e # Exit on error

${CONDA_PREFIX}/bin/R --vanilla -q -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR.results@f92c22a", dependencies = TRUE, upgrade = "never");'

${CONDA_PREFIX}/bin/R --vanilla -q -e 'options(timeout=9999999); devtools::install_github("mg14/mg14", dependencies = TRUE, upgrade = "never")'

${CONDA_PREFIX}/bin/R --vanilla -q -e 'options(timeout=9999999); devtools::install_github("gerstung-lab/MutationTimeR@e4e266a", dependencies = TRUE, upgrade = "never")'