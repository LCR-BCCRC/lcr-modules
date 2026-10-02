#!/usr/bin/env python3
"""
Extract nucleotide and amino acid FASTA from MiXCR exportClonesPretty output.

For each clone, every assembled target (Target0, Target1, ...) is written as a
separate FASTA record:
  >cloneId_N_readCount_X_readFraction_Y_targetZ

Nucleotide FASTA: TargetZ segment sequences concatenated in position order,
  alignment gaps (-) stripped.
Amino acid FASTA: translated codes above each TargetZ block concatenated in
  block order; partial-codon markers (_) excluded, stop codons (*) retained.

Translation notation follows MiXCR conventions:
  https://mixcr.com/mixcr/reference/ref-translation-rules/
"""

import re
import argparse

# TargetN nucleotide line: "    Target0  0 ACGT... 79   Score (hit score)"
_TARGET_RE = re.compile(
    r'^\s+Target(\d+)\s+(\d+)\s+([A-Za-z\-]+)\s+\d+\s+Score'
)

# AA annotation line: starts with whitespace, contains uppercase single-letter
# codes (A-Z, *, X) separated by spaces; underscore _ marks partial codons.
# Requires at least 2 codes to avoid false positives.
_AA_LINE_RE = re.compile(r'^\s+[A-Z_*X](\s+[A-Z_*X])+\s*$')

_CLONE_HEADER_RE = re.compile(r'^>>> Clone id:\s*(\S+)')
_ABUNDANCE_RE = re.compile(
    r'^>>> Abundance, reads \(fraction\):\s*(\S+)\s+\((\S+)\)'
)


def _extract_aa_codes(line):
    # [A-Z*] captures standard AA codes and stop codons; _ is excluded
    return re.findall(r'[A-Z*]', line)


def parse_pretty_file(filepath):
    with open(filepath) as fh:
        text = fh.read()

    clones = []
    # Split on clone boundaries (lookahead keeps the marker on the right chunk)
    blocks = re.split(r'(?=^>>> Clone id:)', text, flags=re.MULTILINE)

    for block in blocks:
        if not block.strip() or not block.startswith('>>> Clone id:'):
            continue

        clone_id = read_count = read_fraction = None
        for line in block.split('\n'):
            m = _CLONE_HEADER_RE.match(line)
            if m:
                clone_id = m.group(1)
            m = _ABUNDANCE_RE.match(line)
            if m:
                read_count, read_fraction = m.group(1), m.group(2)

        if clone_id is None:
            continue

        targets_nt = {}  # tnum -> [(start, seq), ...]
        targets_aa = {}  # tnum -> [aa_code, ...]
        aa_buffer = []   # AA codes accumulated since last Target line

        for line in block.split('\n'):
            if _AA_LINE_RE.match(line):
                aa_buffer.extend(_extract_aa_codes(line))
                continue
            m = _TARGET_RE.match(line)
            if m:
                tnum = int(m.group(1))
                start = int(m.group(2))
                seq = m.group(3).replace('-', '')  # strip alignment gaps
                targets_nt.setdefault(tnum, []).append((start, seq))
                if aa_buffer:
                    targets_aa.setdefault(tnum, []).extend(aa_buffer)
                    aa_buffer = []

        nt_seqs = {
            t: ''.join(s for _, s in sorted(segs))
            for t, segs in targets_nt.items()
        }
        aa_seqs = {
            t: ''.join(codes)
            for t, codes in targets_aa.items()
        }

        clones.append({
            'clone_id': clone_id,
            'read_count': read_count,
            'read_fraction': read_fraction,
            'nt_seqs': nt_seqs,
            'aa_seqs': aa_seqs,
        })

    return clones


def main():
    parser = argparse.ArgumentParser(
        description='Extract NT and AA FASTA from MiXCR exportClonesPretty output'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='exportClonesPretty .txt file')
    parser.add_argument('--output_nt', required=True,
                        help='output nucleotide FASTA')
    parser.add_argument('--output_aa', required=True,
                        help='output amino acid FASTA')
    args = parser.parse_args()

    clones = parse_pretty_file(args.input)

    with open(args.output_nt, 'w') as nt_out, \
         open(args.output_aa, 'w') as aa_out:
        for clone in clones:
            base = (
                f"cloneId_{clone['clone_id']}"
                f"_readCount_{clone['read_count']}"
                f"_readFraction_{clone['read_fraction']}"
            )
            for tnum in sorted(clone['nt_seqs']):
                seq = clone['nt_seqs'][tnum]
                if seq:
                    nt_out.write(f">{base}_target{tnum}\n{seq}\n")
            for tnum in sorted(clone['aa_seqs']):
                seq = clone['aa_seqs'][tnum]
                if seq:
                    aa_out.write(f">{base}_target{tnum}\n{seq}\n")


if __name__ == '__main__':
    main()
