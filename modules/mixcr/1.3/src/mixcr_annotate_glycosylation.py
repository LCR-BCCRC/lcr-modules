#!/usr/bin/env python3
"""
Annotate N-linked glycosylation sites (NxS/T, x != Pro) in MiXCR-assembled
immunoglobulin amino acid sequences using IMGT unique numbering (ANARCI).

Input is the amino acid FASTA produced by mixcr_to_fasta.py (--output_aa).

Because no germline reference is available, glycosylation sites cannot be
definitively classified as SHM-acquired.  Sites located in CDR1, CDR2, or
CDR3 are flagged separately via `potential_ags` as they are more likely to
represent somatic mutations than FWR sites; however, this is not equivalent
to the germline-comparison-based AGS classification performed by
annotate_glycosylation.py in the annotate_bcr module.

MiXCR AA sequences contain lowercase characters for germline-imputed regions
and uppercase for observed sequence.  All characters are uppercased before
ANARCI numbering.  Sequences are truncated at the first stop codon (*).

Requires ANARCI (bioconda: conda install -c bioconda anarci).

Output TSV columns:
  sequence_id                   FASTA header ID
  aa_sequence                   Cleaned AA sequence submitted to ANARCI
  num_glycosylation_sites       Total count of NxS/T motifs (x != Pro)
  glycosylation_imgt_positions  IMGT positions of N residues, comma-separated
  glycosylation_motifs          NxS/T triplets, comma-separated
  glycosylation_imgt_regions    IMGT FWR/CDR region of each site, comma-separated
  num_glycosylation_sites_cdr   Sites located in CDR1/CDR2/CDR3
  potential_ags                 POS/NEG: >=1 CDR glycosylation site present.
                                Suggestive of an acquired glycosylation site but
                                germline comparison is not available to confirm.
"""

import argparse
import csv
import re
import sys

try:
    from anarci import anarci
except ImportError:
    sys.exit("ERROR: anarci is not installed. Install via: conda install -c bioconda anarci")


# IMGT V-domain region boundaries (Lefranc et al.), inclusive
IMGT_REGION_BOUNDARIES = [
    (1,   26,  "FWR1"),
    (27,  38,  "CDR1"),
    (39,  55,  "FWR2"),
    (56,  65,  "CDR2"),
    (66,  104, "FWR3"),
    (105, 117, "CDR3"),
    (118, 128, "FWR4"),
]

FIELDS = [
    "sequence_id",
    "aa_sequence",
    "num_glycosylation_sites",
    "glycosylation_imgt_positions",
    "glycosylation_motifs",
    "glycosylation_imgt_regions",
    "num_glycosylation_sites_cdr",
    "potential_ags",
]


def _imgt_pos_str(pos_tuple):
    num, ins = pos_tuple
    return f"{num}{ins.strip()}" if ins.strip() else str(num)


def _imgt_region(pos_tuple):
    num, _ = pos_tuple
    for lo, hi, name in IMGT_REGION_BOUNDARIES:
        if lo <= num <= hi:
            return name
    return "NA"


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


def prepare_aa(raw_seq):
    """Uppercase, truncate at first stop codon, strip non-IUPAC AA characters."""
    seq = raw_seq.upper()
    stop = seq.find("*")
    if stop != -1:
        seq = seq[:stop]
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWXY]", "", seq)


def number_with_imgt(aa_seq):
    """Run ANARCI IMGT numbering; return [(pos_tuple, aa), ...] or None."""
    try:
        results, _, _ = anarci([("seq", aa_seq)], scheme="imgt", output=False)
    except Exception as exc:
        print(f"WARNING: ANARCI failed (len={len(aa_seq)}): {exc}", file=sys.stderr)
        return None
    if results[0] is None:
        return None
    return [(pos, aa) for pos, aa in results[0][0][0] if aa != "-"]


def find_glycosylation_sites(numbered):
    """Return [(imgt_pos_str, motif, region), ...] for all NxS/T motifs."""
    sites = []
    for i in range(len(numbered) - 2):
        n_pos, n_aa  = numbered[i]
        _,     x_aa  = numbered[i + 1]
        _,     st_aa = numbered[i + 2]
        if n_aa == "N" and x_aa != "P" and st_aa in ("S", "T"):
            sites.append((_imgt_pos_str(n_pos), f"N{x_aa}{st_aa}", _imgt_region(n_pos)))
    return sites


def _sites_str(sites):
    return ",".join(s[0] for s in sites) if sites else "NA"

def _motifs_str(sites):
    return ",".join(s[1] for s in sites) if sites else "NA"

def _regions_str(sites):
    return ",".join(s[2] for s in sites) if sites else "NA"


def main():
    parser = argparse.ArgumentParser(
        description="Annotate glycosylation sites in MiXCR amino acid FASTA"
    )
    parser.add_argument("--fasta",  required=True, help="Input amino acid FASTA")
    parser.add_argument("--output", required=True, help="Output TSV")
    args = parser.parse_args()

    seqs = read_fasta(args.fasta)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id, raw_seq in seqs.items():
            aa_seq = prepare_aa(raw_seq)

            if not aa_seq:
                writer.writerow({
                    "sequence_id": seq_id,
                    "aa_sequence": raw_seq,
                    "num_glycosylation_sites": 0,
                    "glycosylation_imgt_positions": "NA",
                    "glycosylation_motifs": "NA",
                    "glycosylation_imgt_regions": "NA",
                    "num_glycosylation_sites_cdr": 0,
                    "potential_ags": "NEG",
                })
                continue

            numbered  = number_with_imgt(aa_seq)
            sites     = find_glycosylation_sites(numbered) if numbered else []
            cdr_sites = [s for s in sites if s[2].startswith("CDR")]

            writer.writerow({
                "sequence_id": seq_id,
                "aa_sequence": aa_seq,
                "num_glycosylation_sites": len(sites),
                "glycosylation_imgt_positions": _sites_str(sites),
                "glycosylation_motifs": _motifs_str(sites),
                "glycosylation_imgt_regions": _regions_str(sites),
                "num_glycosylation_sites_cdr": len(cdr_sites),
                "potential_ags": "POS" if cdr_sites else "NEG",
            })


if __name__ == "__main__":
    main()
