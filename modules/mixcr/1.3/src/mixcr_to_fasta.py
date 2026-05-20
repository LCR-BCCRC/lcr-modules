#!/usr/bin/env python

# Extracts FR1, CDR1, FR2, CDR2, FR3, CDR3 and FR4 sequences from MiXCR results and joints into fasta format
import os
import argparse
import sys 

# Parse arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',required=True)
    parser.add_argument('-o','--output',required=True)
    parser.add_argument('-s', '--sequence', required=True)
    args = parser.parse_args()

    return args

def main(args):
    mixcr_input = args.input
    fasta_output = args.output
    output_file = os.path.abspath(fasta_output)
    if os.path.exists(output_file):
        sys.exit(f"{fasta_output} already exists.")
    else:
        out = open(fasta_output,'w')
    seq_info_file = args.sequence
    seq_info = open(seq_info_file, 'a')
    seq_header = "cloneId\tcloneFraction\tcloneCount\tnumMissing\tregionsMissing\n"
    seq_info.write(seq_header)
    # Get column names and positions
    positions = {}
    count = 0
    with open(mixcr_input,'r') as handle:
        for line in handle:
            line = line.strip('\n')
            count += 1
            columns = line.split('\t')
            if count == 1:
                pos = 0
                for feature in columns:
                    positions[feature] = pos
                    pos += 1
            # Extract clone Id, clone fraction, clone count
            if count > 1:
                cloneId = columns[positions["cloneId"]]
                cloneFraction = columns[positions["cloneFraction"]]
                cloneCount = columns[positions["cloneCount"]]
                nSeqFR1 = columns[positions["nSeqFR1"]]
                nSeqCDR1 = columns[positions["nSeqCDR1"]]
                nSeqFR2 = columns[positions["nSeqFR2"]]
                nSeqCDR2 = columns[positions["nSeqCDR2"]]
                nSeqFR3 = columns[positions["nSeqFR3"]]
                nSeqCDR3 = columns[positions["nSeqCDR3"]]
                nSeqFR4 = columns[positions["nSeqFR4"]]
                header = f">cloneId_{cloneId}_cloneFraction_{cloneFraction}_cloneCount_{cloneCount}\n"
                out.write(header)
                sequence = f"{nSeqFR1}{nSeqCDR1}{nSeqFR2}{nSeqCDR2}{nSeqFR3}{nSeqCDR3}{nSeqFR4}\n"
                out.write(sequence)
                num_missing = 0
                regions_missing = []
                position = 0
                string_regions = ["FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4"]
                for region in [nSeqFR1,nSeqCDR1,nSeqFR2,nSeqCDR2,nSeqFR3,nSeqCDR3,nSeqFR4]:
                    if region == "":
                        num_missing += 1
                        regions_missing.append(string_regions[position])
                    position += 1
                regions_missing_string = ",".join(regions_missing)
                num_missing_string = str(num_missing)
                seq_info_line = f"{cloneId}\t{cloneFraction}\t{cloneCount}\t{num_missing_string}\t{regions_missing_string}\n"
                seq_info.write(seq_info_line)
    out.close()
    seq_info.close()

if __name__ == "__main__":
    args = parse_args()
    main(args)