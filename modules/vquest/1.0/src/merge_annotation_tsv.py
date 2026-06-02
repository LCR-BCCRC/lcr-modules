#!/usr/bin/env python3
"""
Left-join an annotation TSV onto a base TSV by a shared key column.
Rows in the base table that have no match in the annotation table get 'N/A'
for every annotation column.

When --base_key and --annot_key differ, annotation key values are parsed
through _extract_mixcr_clone_id to handle MiXCR compound FASTA headers
("cloneId_N_readFraction_F_readCount_C") that encode the base join key.
When both key columns are the same no transformation is applied.
"""

import argparse
import csv
import re


def _extract_mixcr_clone_id(sequence_id):
    """Parse cloneId from a MiXCR compound sequence_id.

    mixcr_to_fasta.py writes FASTA headers as:
        cloneId_{N}_readFraction_{F}_readCount_{C}
    Returns the cloneId string, or the original value if the format does not match.
    """
    m = re.match(r"cloneId_(\S+?)_readFraction_", sequence_id)
    return m.group(1) if m else sequence_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base",       required=True, help="base TSV (all rows kept)")
    parser.add_argument("--annotation", required=True, help="annotation TSV to join in")
    parser.add_argument("--base_key",   default="sequence_id", help="join key column in base")
    parser.add_argument("--annot_key",  default="sequence_id", help="join key column in annotation")
    parser.add_argument("--output",     required=True, help="output TSV path")
    args = parser.parse_args()

    # Transform annotation keys when the key columns differ (e.g. base uses cloneId
    # but the annotation encodes it in a compound sequence_id).
    transform_key = args.base_key != args.annot_key

    annot, annot_cols = {}, []
    with open(args.annotation) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        annot_cols = [c for c in reader.fieldnames if c != args.annot_key]
        for row in reader:
            raw_key = row[args.annot_key]
            key = _extract_mixcr_clone_id(raw_key) if transform_key else raw_key
            annot[key] = {c: row[c] for c in annot_cols}

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
