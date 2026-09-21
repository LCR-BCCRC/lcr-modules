#!/usr/bin/env python3
"""
Generate IgSeqR per-chain reports from kallisto abundance and extracted transcripts.

Outputs:
  report.tsv          - all transcripts with TPM and sequence, sorted by TPM descending
  dominant_report.tsv - top N transcripts by TPM
  tpm_fasta           - FASTA of top N transcripts
"""

import argparse
import csv
import sys


def _read_fasta(path):
    seqs = {}
    current_id, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(parts)
    return seqs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--abundance",       required=True, help="kallisto abundance.tsv")
    parser.add_argument("--transcripts",     required=True, help="per-chain transcripts FASTA")
    parser.add_argument("--sample_id",       required=True)
    parser.add_argument("--report",          required=True)
    parser.add_argument("--dominant_report", required=True)
    parser.add_argument("--tpm_fasta",       required=True)
    parser.add_argument("--n_dominant",      type=int, default=5,
                        help="number of top-TPM transcripts to retain (default: 5)")
    args = parser.parse_args()

    seqs = _read_fasta(args.transcripts)

    rows = []
    with open(args.abundance) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid = row["target_id"]
            rows.append({
                "sample_id":     args.sample_id,
                "transcript_id": tid,
                "length":        row["length"],
                "eff_length":    row["eff_length"],
                "est_counts":    row["est_counts"],
                "tpm":           row["tpm"],
                "sequence":      seqs.get(tid, ""),
            })

    rows.sort(key=lambda r: float(r["tpm"]), reverse=True)

    header = ["sample_id", "transcript_id", "length", "eff_length", "est_counts", "tpm", "sequence"]

    with open(args.report, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    dominant = rows[: args.n_dominant]

    with open(args.dominant_report, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(dominant)

    with open(args.tpm_fasta, "w") as fh:
        for r in dominant:
            if r["sequence"]:
                fh.write(f">{args.sample_id}_{r['transcript_id']}\n{r['sequence']}\n")


if __name__ == "__main__":
    main()
