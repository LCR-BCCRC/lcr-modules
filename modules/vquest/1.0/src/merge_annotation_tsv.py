#!/usr/bin/env python3
"""
Left-join an annotation TSV onto a base TSV by a shared key column.
Rows in the base table that have no match in the annotation table get 'N/A'
for every annotation column.
"""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base",       required=True, help="base TSV (all rows kept)")
    parser.add_argument("--annotation", required=True, help="annotation TSV to join in")
    parser.add_argument("--base_key",   default="sequence_id", help="join key column in base")
    parser.add_argument("--annot_key",  default="sequence_id", help="join key column in annotation")
    parser.add_argument("--output",     required=True, help="output TSV path")
    args = parser.parse_args()

    annot, annot_cols = {}, []
    with open(args.annotation) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        annot_cols = [c for c in reader.fieldnames if c != args.annot_key]
        for row in reader:
            annot[row[args.annot_key]] = {c: row[c] for c in annot_cols}

    with open(args.base) as fh_in, open(args.output, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")
        header = list(reader.fieldnames) + annot_cols
        writer = csv.DictWriter(fh_out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            row.update(annot.get(row[args.base_key], {c: "N/A" for c in annot_cols}))
            writer.writerow(row)


if __name__ == "__main__":
    main()
