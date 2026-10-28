#!/usr/bin/env python

import argparse
import os
import pysam
from collections import Counter
from sortedcontainers import SortedDict


class geom_coord():
    # This is REALLY stupid, but because genomic coordinates are generally a mix of strings and integers
    # ex. chr2:145320
    # we need to set a custom sort order to ensure the output files are sorted correctly
    def __init__(self, chrom, pos):
        self.coord = chrom + ":" + str(pos)
        self.chrom = chrom
        self.pos = pos

    def __hash__(self):
        return hash(self.coord)

    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos

    def __lt__(self, other):
        # Sort by chromosome first, then pos
        if self.chrom == other.chrom:
            return self.pos < other.pos
        return self.chrom < other.chrom


def is_valid_file(filepath, parser):
    if not os.path.exists(filepath):
        raise parser.error(f"Unable to locate \'{filepath}\': No such file or directory")
    else:
        return filepath


def is_valid_fraction(fraction, parser):
    fraction = float(fraction)
    if 0 <= fraction <= 1:
        return fraction
    else:
        raise parser.error(f"Invalid -f/--fraction {fraction} provided. Must be between 0 and 1")

def get_args():
    parser = argparse.ArgumentParser(description="Generates a TSV file of all genomic coordinates with a high proportion of masked bases")
    parser.add_argument("-i", "--input", metavar="BAM/CRAM", required=True, type=lambda x: is_valid_file(x, parser), help="Input BAM/CRAM file. Must be coordinate sorted and indexed")
    parser.add_argument("-r", "--regions", metavar="BED", default=None, type=lambda x: is_valid_file(x, parser), help="A BED file containing regions to restrict to")
    parser.add_argument("-o", "--output", metavar="TSV", required=True, type=str, help="Output TSV file listing all positions with greater than -f/--fraction of N bases")
    parser.add_argument("-f", "--fraction", metavar="FLOAT", default=0.1, type=lambda x: is_valid_fraction(x, parser), help="Flag positions with above this fraction of N bases")
    parser.add_argument("-c", "--count", metavar="INT", default=8, type=int, help="Require at least this many N's to mask a position")

    return parser.parse_args()


def load_regions(bedfile):
    """
    Load a set of regions from a provided BED file

    This function performs additional checks to ensure the BED entries are valid
    These include:
       - Complete line (BED3 minimum)
       - Valid coordinates (i.e. column 2 and 3 are actually positive numbers)
       - Start < End
       - Remove duplicate entries

    Returns a list of the regions in a <chr>:<start>-<end> format

    """
    coords = list()
    total_size = 0
    with open(bedfile) as f:
        i = 0  # Line counter
        for line in f:
            i += 1
            # Skip comment lines, if they exist
            if line.startswith("#"):
                continue

            line = line.rstrip("\n").rstrip("\r")
            cols = line.split("\t")
            try:
                chrom = cols[0]
                start = int(cols[1])
                end = int(cols[2])
                length = end - start
                # Sanity check coordinates
                if start < 0:
                    raise AttributeError(f"Start coordinate {start} is less than 0")
                if end <= start:
                    raise AttributeError(f"Start coordinate {start} is after end coordinate {end}")
            except (TypeError, AttributeError) as e:  # Start and/or end coordinates are not numbers!
                raise AttributeError(f"Problem encountered when processing line {i} of \'{bedfile}\': \'{line}\'") from e
            except IndexError as e:  # Truncated line
                raise AttributeError(f"Unable to process line {i} of \'{bedfile}\' as it appears to be truncated: \'{line}\'") from e

            # Coordinates are valid. Store them as a region for the pileup
            coord = chrom + ":" + str(start) + "-" + str(end)
            coords.append(coord)
            total_size += length

    # Remove duplicates
    return list(set(coords)), str(total_size)


def main(args=None):
    if args is None:
        args = get_args()

    # Open input file
    infile = pysam.AlignmentFile(args.input, threads=3)

    # Load regions, if they are provided
    if args.regions is not None:
        regions, region_size = load_regions(args.regions)
    else:  # No regions were provided. Process all regions with coverage
        regions = [None]
        region_size = "All_covered"

    # Perform the pileup
    n_fracs = SortedDict()
    i = 0  # Counter for status messages
    for region in regions:
        pileup = infile.pileup(region=region, truncate=True, max_depth=100000, stepper="nofilter", ignore_overlaps=False, flag_filter=-1, ignore_orphans=False, min_base_quality=0, compute_baq=False)
        for position in pileup:
            i += 1
            coord_key = geom_coord(position.reference_name, position.reference_pos)
            # Get the bases at this position
            bases = position.get_query_sequences()
            # Count the number of each bases
            base_count = Counter(list(x.upper() for x in bases))

            # Print status update message
            if i % 10000 == 0:
                print(f"Processing position {i}/{region_size}: {coord_key.coord}")
            # Get total depth and the number of Ns
            depth = sum(x for x in base_count.values())
            try:
                n_depth = base_count["N"]
            except KeyError:
                # No N's at this position
                n_depth = 0
            try:
                frac_n = n_depth / depth
            except ZeroDivisionError:
                # This shouldn't happen, as there will be no bases to pileup at this position?
                # Might happen if the htslib API changes at some point to auto-drop/ignore certain reads
                # For safety
                frac_n = 0

            # Using a dictionary, as this will drop duplicate positions via overwritting them
            # The overwritten position should have the same stats
            # This can occur if overlapping regions are provided
            n_fracs[coord_key] = [position.reference_name, str(position.reference_pos), str(position.reference_pos + 1), str(n_depth), str(depth), str(frac_n)]

    # Now that we have performed the pileup, filter for positions with too many N's and output them
    with open(args.output, "w") as o:
        for attributes in n_fracs.values():
            # Check fraction of N's
            if float(attributes[5]) >= args.fraction and int(attributes[3]) > args.count:
                o.write("\t".join(attributes))
                o.write(os.linesep)


if __name__ == "__main__":
    main()
