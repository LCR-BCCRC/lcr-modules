#!env bash
set -o pipefail
set -e # Exit on error

#get path to broken file and overwrite it with one we know works

broken_bin=`which faToTwoBit`

#replace broken file with the one that works
curl https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit > $broken_bin

