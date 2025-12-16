#!env bash
set -o pipefail
set -e # Exit on error

R --vanilla -q -e 'devtools::install_github("VanLoo-lab/ascat/ASCAT", dependencies = TRUE, upgrade = "never")'
R --vanilla -q -e 'devtools::install_github("Wedge-Oxford/battenberg", dependencies = FALSE, upgrade = "never")'