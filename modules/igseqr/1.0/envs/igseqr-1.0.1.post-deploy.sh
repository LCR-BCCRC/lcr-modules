#!env bash
set -o pipefail
set -e # Exit on error

# Clone IgSeqR from GitHub and run its setup.sh, which installs the igseqr
# script into $CONDA_PREFIX/bin (the active conda environment).
TMPDIR=$(mktemp -d)
git clone https://github.com/ForconiLab/IgSeqR.git "$TMPDIR/IgSeqR"
cd "$TMPDIR/IgSeqR"
bash setup.sh
cd -
rm -rf "$TMPDIR"
