#!/usr/bin/env python3
"""
Extract FR1-FR4 nucleotide (and optionally amino acid) FASTA from a MiXCR
clonotype TSV.

Column detection is automatic:
  - NT:  prefers nSeqImputed{FR1..FR4} (produced with --impute-germline-on-export),
         falls back to nSeq{FR1..FR4}.
  - AA:  prefers aaSeqImputed{FR1..FR4}, falls back to aaSeq{FR1..FR4}.
         If neither is present, --output_aa is left empty.
"""

import argparse
import csv
import sys

REGIONS = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]


def _detect_prefix(fieldnames, imputed, plain):
    """Return the first prefix for which <prefix>FR1 exists in fieldnames."""
    for p in (imputed, plain):
        if f"{p}FR1" in fieldnames:
            return p
    return None


def _clean_nt(val):
    return "" if val in ("region_not_covered", "") else val


def _clean_aa(val):
    """Strip partial-codon markers (_) and discard region_not_covered."""
    return "" if val in ("region_not_covered", "") else val.replace("_", "")


def main():
    parser = argparse.ArgumentParser(
        description="Convert MiXCR clonotype TSV to FASTA format"
    )
    parser.add_argument("-i", "--input",     required=True, help="MiXCR clonotype TSV")
    parser.add_argument("-o", "--output",    required=True, help="output nucleotide FASTA")
    parser.add_argument("-s", "--sequence",  required=True, help="output seq-info TSV")
    parser.add_argument("-a", "--output_aa", default=None,  help="output amino acid FASTA")
    args = parser.parse_args()

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    nt_prefix = _detect_prefix(fieldnames, "nSeqImputed", "nSeq")
    aa_prefix = _detect_prefix(fieldnames, "aaSeqImputed", "aaSeq")

    if nt_prefix is None:
        sys.exit(
            "ERROR: could not find nSeqImputedFR1 or nSeqFR1 columns in "
            f"{args.input}. Check that the MiXCR export includes per-region "
            "nucleotide sequences."
        )

    nt_cols = [f"{nt_prefix}{r}" for r in REGIONS]
    aa_cols = [f"{aa_prefix}{r}" for r in REGIONS] if aa_prefix else []

    with open(args.output, "w") as nt_out, \
         open(args.sequence, "w") as si_out:

        aa_out = open(args.output_aa, "w") if args.output_aa else None

        si_out.write("cloneId\treadFraction\treadCount\tnumMissing\tregionsMissing\n")

        for row in rows:
            clone_id      = row["cloneId"]
            read_fraction = row["readFraction"]
            read_count    = row["readCount"]
            header        = f">cloneId_{clone_id}_readFraction_{read_fraction}_readCount_{read_count}\n"

            nt_seqs = [_clean_nt(row.get(c, "")) for c in nt_cols]
            nt_out.write(header)
            nt_out.write("".join(nt_seqs) + "\n")

            missing = [r for r, s in zip(REGIONS, nt_seqs) if s == ""]
            si_out.write(
                f"{clone_id}\t{read_fraction}\t{read_count}\t"
                f"{len(missing)}\t{','.join(missing)}\n"
            )

            if aa_out and aa_cols:
                aa_seqs = [_clean_aa(row.get(c, "")) for c in aa_cols]
                aa_seq  = "".join(aa_seqs)
                if aa_seq:
                    aa_out.write(header)
                    aa_out.write(aa_seq + "\n")

        if aa_out:
            aa_out.close()


if __name__ == "__main__":
    main()
