#!/usr/bin/env python3
"""
Subset a MiXCR clonotype TSV to the top N clones by readFraction and write
a matching FASTA and seq_info file. MiXCR outputs clones in descending
readFraction order, so taking the first N rows gives the most abundant clones.
"""

import argparse
import csv


REGIONS = ["nSeqFR1", "nSeqCDR1", "nSeqFR2", "nSeqCDR2", "nSeqFR3", "nSeqCDR3", "nSeqFR4"]
REGION_LABELS = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",        required=True, help="MiXCR clonotype TSV")
    parser.add_argument("-o", "--output_tsv",   required=True, help="output dominant TSV")
    parser.add_argument("-f", "--output_fasta",  required=True, help="output dominant FASTA")
    parser.add_argument("-s", "--output_seq_info", required=True, help="output dominant seq_info")
    parser.add_argument("-n", "--n_dominant",   type=int, default=5,
                        help="number of top clones to retain (default: 5)")
    args = parser.parse_args()

    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    dominant = rows[: args.n_dominant]

    with open(args.output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(dominant)

    with open(args.output_fasta, "w") as fa, open(args.output_seq_info, "w") as si:
        si.write("cloneId\treadFraction\treadCount\tnumMissing\tregionsMissing\n")
        for row in dominant:
            clone_id      = row["cloneId"]
            read_fraction = row["readFraction"]
            read_count    = row["readCount"]

            seqs = [
                "" if row.get(r, "") == "region_not_covered" else row.get(r, "")
                for r in REGIONS
            ]
            seq = "".join(seqs)
            fa.write(f">cloneId_{clone_id}_readFraction_{read_fraction}_readCount_{read_count}\n{seq}\n")

            missing = [label for label, s in zip(REGION_LABELS, seqs) if s == ""]
            si.write(f"{clone_id}\t{read_fraction}\t{read_count}\t{len(missing)}\t{','.join(missing)}\n")


if __name__ == "__main__":
    main()
