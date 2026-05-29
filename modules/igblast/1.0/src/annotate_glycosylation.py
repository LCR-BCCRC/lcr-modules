#!/usr/bin/env python3
"""
Annotate acquired N-linked glycosylation sites (NxS/T, x != Pro) in IG/TCR
sequences with IMGT unique numbering (Lefranc et al.).

Two input modes (auto-detected from --source_tsv):

  igseqr / default
    Translates each cDNA sequence in all six reading frames (3 forward +
    3 reverse complement) and selects the longest stop-codon-free ORF.

  MiXCR (triggered when --source_tsv is provided and "mixcr" is in the path)
    Assembles the amino acid sequence from aaSeqFR1..aaSeqFR4 columns in the
    clonotype TSV, skipping fields whose value is "region_not_covered".
    The FASTA is still required to supply sequence IDs; each FASTA entry's
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


# ── Translation ────────────────────────────────────────────────────────────────

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _translate_frame(nt_seq, frame):
    aa = []
    pos = frame
    while pos + 3 <= len(nt_seq):
        aa.append(CODON_TABLE.get(nt_seq[pos:pos + 3].upper(), "X"))
        pos += 3
    return "".join(aa)


def best_translation(nt_seq):
    """Return the longest stop-codon-free segment across all six reading frames."""
    rc = nt_seq.translate(_COMPLEMENT)[::-1]
    best = ""
    for seq in (nt_seq, rc):
        for frame in range(3):
            for segment in _translate_frame(seq, frame).split("*"):
                if len(segment) > len(best):
                    best = segment
    return best


# ── MiXCR AA assembly ──────────────────────────────────────────────────────────

_AA_COLS = [
    "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2",
    "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4",
]

_MIXCR_CLONE_RE = re.compile(r"cloneId_(\S+?)_readFraction_")


def _parse_clone_id(sequence_id):
    """Extract the bare cloneId from a MiXCR FASTA header."""
    m = _MIXCR_CLONE_RE.match(sequence_id)
    return m.group(1) if m else None


def load_mixcr_aa(tsv_path):
    """
    Build {cloneId: assembled_aa_sequence} from a MiXCR clonotype TSV.
    Columns with value "region_not_covered" are skipped.
    """
    aa_by_clone = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            segments = [
                row[col] for col in _AA_COLS
                if col in row and row[col] and row[col] != "region_not_covered"
            ]
            aa_seq = "".join(segments)
            if aa_seq:
                aa_by_clone[row["cloneId"]] = aa_seq
    return aa_by_clone


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
                        help="Input cDNA FASTA. Used as AA source for igseqr "
                             "(6-frame translation) and as sequence-ID source "
                             "for MiXCR (AA comes from --source_tsv instead).")
    parser.add_argument("--source_tsv", default=None,
                        help="MiXCR clonotype TSV with aaSeq* columns. When "
                             "provided and 'mixcr' is in the path, AA sequences "
                             "are assembled from the TSV rather than translated.")
    parser.add_argument("--output",     required=True, help="Output TSV path")
    args = parser.parse_args()

    fasta_seqs = read_fasta(args.fasta)

    mixcr_mode = bool(args.source_tsv and "mixcr" in args.source_tsv.lower())
    mixcr_aa = load_mixcr_aa(args.source_tsv) if mixcr_mode else {}

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()

        for seq_id, nt_seq in fasta_seqs.items():

            # ── Resolve amino acid sequence ────────────────────────────────
            if mixcr_mode:
                clone_id = _parse_clone_id(seq_id)
                aa_seq = mixcr_aa.get(clone_id, "") if clone_id else ""
                if not aa_seq:
                    aa_seq = best_translation(nt_seq)  # fallback
            else:
                aa_seq = best_translation(nt_seq)

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