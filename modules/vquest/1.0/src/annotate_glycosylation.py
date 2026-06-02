#!/usr/bin/env python3
"""
Annotate acquired N-linked glycosylation sites (NxS/T, x != Pro) in IG/TCR
sequences with IMGT unique numbering (Lefranc et al.).

Reads sequence_alignment_aa from the V-QUEST AIRR TSV (--source_tsv), keyed
by sequence_id. Works for both igseqr and MiXCR input sources.

IMGT numbering is assigned by ANARCI (bioconda: conda install -c bioconda anarci).

Output TSV columns:
  sequence_id                  FASTA header ID
  aa_sequence                  Amino acid sequence used for numbering
  num_glycosylation_sites      Count of NxS/T motifs (x != Pro)
  glycosylation_imgt_positions IMGT positions of the N residues, comma-separated
  glycosylation_motifs         Corresponding NxS/T triplets, comma-separated
"""

import argparse
import csv
import sys

try:
    from anarci import anarci
except ImportError:
    sys.exit("ERROR: anarci is not installed. Install via: conda install -c bioconda anarci")


# ── TSV loading ───────────────────────────────────────────────────────────────

def load_aa_by_sequence_id(tsv_path):
    """Build {sequence_id: sequence_alignment_aa} from a V-QUEST AIRR TSV."""
    aa_by_id = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            aa_seq = row.get("sequence_alignment_aa", "").strip()
            if aa_seq:
                aa_by_id[row["sequence_id"]] = aa_seq
    return aa_by_id


# ── IMGT numbering ────────────────────────────────────────────────────────────

def _imgt_pos_str(pos_tuple):
    """Format ANARCI position (number, insertion_letter) as '27' or '111a'."""
    num, ins = pos_tuple
    return f"{num}{ins.strip()}" if ins.strip() else str(num)


def number_with_imgt(aa_seq):
    """
    Align aa_seq to IG/TCR germline HMMs with ANARCI using the IMGT scheme.
    Returns [(pos_tuple, aa), ...] with gap positions removed, or None on failure.
    """
    results, _, _ = anarci([("seq", aa_seq)], scheme="imgt", output=False)
    if results[0] is None:
        return None
    return [(pos, aa) for pos, aa in results[0][0][0] if aa != "-"]


def find_glycosylation_sites(numbered):
    """
    Scan IMGT-numbered residues for acquired N-linked glycosylation motifs.

    Motif: NxS/T where x is any residue except Proline. Checks consecutive
    residues in the ANARCI-numbered list (gaps already removed), so the two
    neighbours are the next and second-next residues in the primary sequence.

    Returns [(imgt_position_string, motif_string), ...] for each site.
    """
    sites = []
    for i in range(len(numbered) - 2):
        _, n_aa  = numbered[i]
        _, x_aa  = numbered[i + 1]
        _, st_aa = numbered[i + 2]
        if n_aa == "N" and x_aa != "P" and st_aa in ("S", "T"):
            sites.append((_imgt_pos_str(numbered[i][0]), f"N{x_aa}{st_aa}"))
    return sites


# ── I/O ───────────────────────────────────────────────────────────────────────

def read_fasta(path):
    seqs = {}
    cur_id, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(parts)
                cur_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(parts)
    return seqs


FIELDS = [
    "sequence_id", "aa_sequence", "num_glycosylation_sites",
    "glycosylation_imgt_positions", "glycosylation_motifs",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",      required=True,
                        help="Input FASTA. Used as the source of sequence IDs.")
    parser.add_argument("--source_tsv", required=True,
                        help="V-QUEST AIRR TSV with a sequence_alignment_aa column.")
    parser.add_argument("--output",     required=True, help="Output TSV path")
    args = parser.parse_args()

    fasta_seqs = read_fasta(args.fasta)
    aa_lookup  = load_aa_by_sequence_id(args.source_tsv)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id in fasta_seqs:
            aa_seq = aa_lookup.get(seq_id, "")

            if not aa_seq:
                writer.writerow({
                    "sequence_id": seq_id, "aa_sequence": "N/A",
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "N/A",
                    "glycosylation_motifs": "N/A",
                })
                continue

            numbered = number_with_imgt(aa_seq)
            if numbered is None:
                writer.writerow({
                    "sequence_id": seq_id, "aa_sequence": aa_seq,
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "N/A",
                    "glycosylation_motifs": "N/A",
                })
                continue

            sites = find_glycosylation_sites(numbered)
            writer.writerow({
                "sequence_id": seq_id,
                "aa_sequence": aa_seq,
                "num_glycosylation_sites": len(sites),
                "glycosylation_imgt_positions": ",".join(s[0] for s in sites) if sites else "N/A",
                "glycosylation_motifs": ",".join(s[1] for s in sites) if sites else "N/A",
            })


if __name__ == "__main__":
    main()
