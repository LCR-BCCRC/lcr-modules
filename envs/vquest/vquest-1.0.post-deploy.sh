#!/usr/bin/env bash
set -euo pipefail

# Install vquest from the lkhilton fork, which fixes compatibility with
# V-QUEST >=3.8.0 (moleculeType parameter and updated response handling).
pip install https://github.com/shawhahnlab/vquest.git@0c01f498ba01f4f289e58155c198f16dc7a277ac