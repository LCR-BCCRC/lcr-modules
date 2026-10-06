#!env bash
set -o pipefail
set -e # Exit on error

# Clone IgSeqR and manually replicate the non-interactive steps from setup.sh:
#   - copy bin/igseqr.sh into $CONDA_PREFIX/bin and symlink as 'igseqr'
#   - copy the bundled data/ directory (IG reference FASTAs) into place
# The large HISAT2 genome index is NOT installed here; it must be provided
# separately via the hisat_ref config input.
TMPDIR=$(mktemp -d)
git clone --branch v1.0.1 --depth 1 https://github.com/ForconiLab/IgSeqR.git "$TMPDIR/IgSeqR"
mkdir -p "$CONDA_PREFIX/bin/data/igseqr"
cp "$TMPDIR/IgSeqR/bin/igseqr.sh" "$CONDA_PREFIX/bin/"
chmod +x "$CONDA_PREFIX/bin/igseqr.sh"
ln -sf "$CONDA_PREFIX/bin/igseqr.sh" "$CONDA_PREFIX/bin/igseqr"
cp -r "$TMPDIR/IgSeqR/data/"* "$CONDA_PREFIX/bin/data/igseqr/"
rm -rf "$TMPDIR"
