#!/bin/bash
# 
# Usage: ./prepare_reference_files.sh /path/to/reference_directory/ [num_cores]

set -euf -o pipefail

if [ -z ${1+x} ]; then 
    echo "Error: First argument must be the /path/to/reference_directory/"
    exit 1
fi

REF_DIR="${1}"
NUM_CORES="${2:-24}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

snakemake --cores "${NUM_CORES}" --use-conda --directory "${SCRIPT_DIR}" \
    --snakefile "${SCRIPT_DIR}/prepare_reference_files.smk" \
    --config reference_directory="${REF_DIR}"
    
