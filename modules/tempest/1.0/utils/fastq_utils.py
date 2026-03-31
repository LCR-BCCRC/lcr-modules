import os
import gzip
import logging
import math
import datetime

def generate_read_group(fastq, sample, config):
    # Parses flowcell, lane, and barcode information from FASTQ read names
    # Uses this information (and config file info) to generate a read group line
    if is_gzipped(fastq):
        readname = gzip.open(fastq, "rt").readline()
    else:
        readname = open(fastq, "r").readline()

    # Parse out the attributes for the read group from the read name
    readname = readname.rstrip("\n").rstrip("\r")
    cols = readname.split(":")

    sequencer = cols[0]
    sequencer = sequencer.replace("@","")
    flowcell = cols[2]
    flowcell = flowcell.split("-")[-1]
    lane = cols[3]
    barcode = cols[-1]

    # From the config (generic and should be consistent between runs)
    description = config["lcr-modules"]["cfDNA_umi_workflow"]["readgroup"]["description"]
    centre = config["lcr-modules"]["cfDNA_umi_workflow"]["readgroup"]["centre"]
    platform = config["lcr-modules"]["cfDNA_umi_workflow"]["readgroup"]["platformunit"]
    platformmodel = config["lcr-modules"]["cfDNA_umi_workflow"]["readgroup"]["platformmodel"]

    # Get the date this sample was generated.
    # This isn't perfect, but it should be relatively close to the sequencing date.
    origpath = os.path.realpath(fastq)
    date = datetime.date.fromtimestamp(os.path.getctime(origpath)).isoformat()  # Get creation date of FASTQ file
    platformunit = flowcell + "-" + lane + ":" + barcode

    readgroup = f"@RG\\tID:{sample}\\tBC:{barcode}\\tCN:{centre}\\tDS:\'{description}\'\\tDT:{date}\\tLB:{sample}\\tPL:{platform}\\tPM:{platformmodel}\\tPU:{platformunit}\\tSM:{sample}"
    return readgroup

def is_gzipped(filepath):
    with open(filepath, "rb") as f:
        magicnum = f.read(2)
        return magicnum == b'\x1f\x8b'  # Magic number of gzipped files

# Adapted from the multiQC module
# Used to calculate the estimated total number of molecules in a library
def estimateLibrarySize(readPairs, uniqueReadPairs):
    """
    Picard calculation to estimate library size
    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L153-L164
    Note: Optical duplicates are contained in duplicates and therefore do not enter the calculation here.
    See also the computation of READ_PAIR_NOT_OPTICAL_DUPLICATES.
     * Estimates the size of a library based on the number of paired end molecules observed
     * and the number of unique pairs observed.
     * <p>
     * Based on the Lander-Waterman equation that states:
     * C/X = 1 - exp( -N/X )
     * where
     * X = number of distinct molecules in library
     * N = number of read pairs
     * C = number of distinct fragments observed in read pairs
    """

    readPairDuplicates = readPairs - uniqueReadPairs

    if readPairs > 0 and readPairDuplicates > 0:

        m = 1.0
        M = 100.0

        if uniqueReadPairs >= readPairs or f(m * uniqueReadPairs, uniqueReadPairs, readPairs) < 0:
            logging.warning("Picard recalculation of ESTIMATED_LIBRARY_SIZE skipped - metrics look wrong")
            return None

        # find value of M, large enough to act as other side for bisection method
        while f(M * uniqueReadPairs, uniqueReadPairs, readPairs) > 0:
            M *= 10.0

        # use bisection method (no more than 40 times) to find solution
        for i in range(40):
            r = (m + M) / 2.0
            u = f(r * uniqueReadPairs, uniqueReadPairs, readPairs)
            if u == 0:
                break
            elif u > 0:
                m = r
            elif u < 0:
                M = r

        return uniqueReadPairs * (m + M) / 2.0
    else:
        return None

def f(x, c, n):
    """
    Picard calculation used when estimating library size
    Taken & translated from the Picard codebase:
    https://github.com/broadinstitute/picard/blob/78ea24466d9bcab93db89f22e6a6bf64d0ad7782/src/main/java/picard/sam/DuplicationMetrics.java#L172-L177
    * Method that is used in the computation of estimated library size.
    """
    return c / x - 1 + math.exp(-n / x)