#!/usr/bin/env python3
"""
Annotate acquired N-linked glycosylation sites (NxS/T, x != Pro) in IG/TCR
sequences with IMGT unique numbering (Lefranc et al.).

Two input modes (auto-detected from --source_tsv):

  igseqr / default
    Reads sequence_alignment_aa from the V-QUEST output TSV, keyed by
    sequence_id.

  MiXCR (triggered when "mixcr" is in the --source_tsv path)
    Reads sequence_alignment_aa from the clonotype TSV, keyed by cloneId.
    The FASTA is required to supply sequence IDs; each FASTA entry's
    cloneId is parsed from its header and looked up in the TSV.

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
import re
import sys

try:
    from anarci import anarci
except ImportError:
    sys.exit("ERROR: anarci is not installed. Install via: conda install -c bioconda anarci")


# ── TSV loading ───────────────────────────────────────────────────────────────

_MIXCR_CLONE_RE = re.compile(r"cloneId_(\S+?)_readFraction_")


def _parse_clone_id(sequence_id):
    """Extract the bare cloneId from a MiXCR FASTA header."""
    m = _MIXCR_CLONE_RE.match(sequence_id)
    return m.group(1) if m else None


def load_aa_by_clone(tsv_path):
    """Build {cloneId: sequence_alignment_aa} from a MiXCR clonotype TSV."""
    aa_by_clone = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            aa_seq = row.get("sequence_alignment_aa", "").strip()
            if aa_seq:
                aa_by_clone[row["cloneId"]] = aa_seq
    return aa_by_clone


def load_aa_by_sequence_id(tsv_path):
    """Build {sequence_id: sequence_alignment_aa} from a V-QUEST output TSV."""
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
                        help="TSV with a sequence_alignment_aa column. For MiXCR "
                             "('mixcr' in path), cloneId is used as the key; "
                             "otherwise sequence_id is used.")
    parser.add_argument("--output",     required=True, help="Output TSV path")
    args = parser.parse_args()

    fasta_seqs = read_fasta(args.fasta)

    mixcr_mode = "mixcr" in args.source_tsv.lower()
    if mixcr_mode:
        aa_lookup = load_aa_by_clone(args.source_tsv)
    else:
        aa_lookup = load_aa_by_sequence_id(args.source_tsv)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id in fasta_seqs:

            # ── Resolve amino acid sequence ────────────────────────────────
            if mixcr_mode:
                clone_id = _parse_clone_id(seq_id)
                aa_seq = aa_lookup.get(clone_id, "") if clone_id else ""
            else:
                aa_seq = aa_lookup.get(seq_id, "")

            if not aa_seq:
                writer.writerow({
                    "sequence_id": seq_id, "aa_sequence": "N/A",
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "N/A",
                    "glycosylation_motifs": "N/A",
                })
                continue

            # ── IMGT numbering ─────────────────────────────────────────────
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
