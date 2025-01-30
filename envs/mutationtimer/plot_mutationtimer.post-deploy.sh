#!env bash
set -o pipefail
set -e # Exit on error

${CONDA_PREFIX}/bin/R --vanilla -q -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMLR.helpers@v1.2", dependencies = TRUE, upgrade = "never")' > /home/sgillis/plot_postdeploy.out
