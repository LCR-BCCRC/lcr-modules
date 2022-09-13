#!/usr/bin/env python


##### ATTRIBUTION #####


# Original Author:  Lakshay Sethi

### Battenberg refrence file corrector ###
# Replaces the placeholder value in the impute_info.txt with the correct path
# where the reference files downloaded are stored.

#
# Usage:
#   python <path/to/src_dir>/reference_correction.py <genome_wildcard>
#
# Notes:
#   This script is intended for use with the Battenberg-1.2 module in LCR-modules.
#   It expects to find the genome build at the input path, following
#   the pattern reference_correction.py {genome_build}. These files should be in the
#   00-inputs subdirectory of the battenberg-1.2 directory present in the results directory.
#
#   The file is made to be present in the src sub directory of the module.
#
#   The sample table should adhere to LCR-modules guidelines.

import os
import sys

cwd = sys.argv[2]

fileIN = open(
    cwd
    + "/impute_info.txt",
    "r",
)
filedata = fileIN.read()
fileIN.close()

newdata = filedata.replace(
    "<path_to_impute_reference_files>",
    cwd
    + "/battenberg_impute_v3",
)

fileOut = open(
    cwd
    + "/impute_info.txt",
    "w",
)
fileOut.write(newdata)
fileOut.close()

