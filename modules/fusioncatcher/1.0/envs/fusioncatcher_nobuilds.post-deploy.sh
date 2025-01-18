#!env bash
set -o pipefail
set -e # Exit on error

#get path to broken file and overwrite it with one we know works

broken_bin=`which faToTwoBit`
echo $broken_bin

# Warning to would-be testers: 
# don't uncomment/test this until the environment is not in use because it will theoretically clobber a file that may be in use!
# curl https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit > $broken_bin

