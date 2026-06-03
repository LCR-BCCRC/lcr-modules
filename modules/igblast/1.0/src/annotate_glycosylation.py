#!/usr/bin/env python3
"""
Annotate acquired N-linked glycosylation sites (NxS/T, x != Pro) in IG/TCR
sequences with IMGT unique numbering (Lefranc et al.).

Reads sequence_alignment_aa and germline_alignment_aa from the igblastn AIRR
TSV (--source_tsv), keyed by sequence_id. Works for both igseqr and MiXCR
input sources.

IMGT numbering is assigned by ANARCI (bioconda: conda install -c bioconda anarci).

Output TSV columns:
  sequence_id                       FASTA header ID
  aa_sequence                       Query amino acid sequence used for numbering
  num_glycosylation_sites           Count of NxS/T motifs in query (x != Pro)
  glycosylation_imgt_positions      IMGT positions of N residues in query, comma-separated
  glycosylation_motifs              Corresponding NxS/T triplets in query, comma-separated
  num_acquired_glycosylation_sites  Sites in query absent from germline (SHM-acquired)
  acquired_glycosylation_imgt_positions  IMGT positions of acquired sites, comma-separated
  germline_aa_sequence              Germline amino acid sequence used for numbering
  num_germline_glycosylation_sites  Count of NxS/T motifs in germline
  germline_glycosylation_imgt_positions  IMGT positions of N residues in germline
  germline_glycosylation_motifs     Corresponding NxS/T triplets in germline
"""

import argparse
import csv
import sys

try:
    from anarci import anarci
except ImportError:
    sys.exit("ERROR: anarci is not installed. Install via: conda install -c bioconda anarci")


# ── TSV loading ───────────────────────────────────────────────────────────────

def load_sequences_by_id(tsv_path):
    """
    Build {sequence_id: {"seq": str, "germline": str}} from an igblastn AIRR TSV.
    Both sequences have alignment gap characters stripped.
    """
    data = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq_id = row["sequence_id"]
            aa_seq  = row.get("sequence_alignment_aa",  "").strip().replace("-", "")
            gl_seq  = row.get("germline_alignment_aa",  "").strip().replace("-", "")
            if aa_seq:
                data[seq_id] = {"seq": aa_seq, "germline": gl_seq}
    return data


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
    Scan IMGT-numbered residues for N-linked glycosylation motifs (NxS/T, x != Pro).

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


def classify_sites(query_sites, germline_sites):
    """
    Return the subset of query_sites whose IMGT position is absent in germline_sites.
    These are SHM-acquired glycosylation sites.
    """
    germline_positions = {pos for pos, _ in germline_sites}
    return [(pos, motif) for pos, motif in query_sites if pos not in germline_positions]


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
    "sequence_id",
    "aa_sequence",
    "num_glycosylation_sites",
    "glycosylation_imgt_positions",
    "glycosylation_motifs",
    "num_acquired_glycosylation_sites",
    "acquired_glycosylation_imgt_positions",
    "germline_aa_sequence",
    "num_germline_glycosylation_sites",
    "germline_glycosylation_imgt_positions",
    "germline_glycosylation_motifs",
]


def _sites_str(sites):
    return ",".join(s[0] for s in sites) if sites else "N/A"

def _motifs_str(sites):
    return ",".join(s[1] for s in sites) if sites else "N/A"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",      required=True,
                        help="Input FASTA. Used as the source of sequence IDs.")
    parser.add_argument("--source_tsv", required=True,
                        help="igblastn AIRR TSV with sequence_alignment_aa and "
                             "germline_alignment_aa columns.")
    parser.add_argument("--output",     required=True, help="Output TSV path")
    args = parser.parse_args()

    fasta_seqs  = read_fasta(args.fasta)
    seq_lookup  = load_sequences_by_id(args.source_tsv)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id in fasta_seqs:
            entry = seq_lookup.get(seq_id)

            if not entry:
                writer.writerow({
                    "sequence_id": seq_id, "aa_sequence": "N/A",
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "N/A",
                    "glycosylation_motifs": "N/A",
                    "num_acquired_glycosylation_sites": 0,
                    "acquired_glycosylation_imgt_positions": "N/A",
                    "germline_aa_sequence": "N/A",
                    "num_germline_glycosylation_sites": 0,
                    "germline_glycosylation_imgt_positions": "N/A",
                    "germline_glycosylation_motifs": "N/A",
                })
                continue

            aa_seq = entry["seq"]
            gl_seq = entry["germline"]

            numbered    = number_with_imgt(aa_seq)
            gl_numbered = number_with_imgt(gl_seq) if gl_seq else None

            query_sites    = find_glycosylation_sites(numbered)    if numbered    else []
            germline_sites = find_glycosylation_sites(gl_numbered) if gl_numbered else []
            acquired_sites = classify_sites(query_sites, germline_sites)

            writer.writerow({
                "sequence_id": seq_id,
                "aa_sequence": aa_seq,
                "num_glycosylation_sites": len(query_sites),
                "glycosylation_imgt_positions": _sites_str(query_sites),
                "glycosylation_motifs": _motifs_str(query_sites),
                "num_acquired_glycosylation_sites": len(acquired_sites),
                "acquired_glycosylation_imgt_positions": _sites_str(acquired_sites),
                "germline_aa_sequence": gl_seq if gl_seq else "N/A",
                "num_germline_glycosylation_sites": len(germline_sites),
                "germline_glycosylation_imgt_positions": _sites_str(germline_sites),
                "germline_glycosylation_motifs": _motifs_str(germline_sites),
            })


if __name__ == "__main__":
    main()
