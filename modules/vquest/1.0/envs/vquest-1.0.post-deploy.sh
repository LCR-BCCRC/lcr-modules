#!/usr/bin/env bash
set -euo pipefail

# Install vquest from the lkhilton fork, which fixes compatibility with
# V-QUEST >=3.8.0 (moleculeType parameter and updated response handling).
pip install git+https://github.com/lkhilton/vquest
