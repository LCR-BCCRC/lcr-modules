#!/usr/bin/env python

# Parses igblastn fmt7 output into a human-readable TSV.
# Compatible with FASTA input from any source (e.g. mixcr 1.3 or igseqr 1.0).

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='igblastn fmt7 output file')
    parser.add_argument('-o', '--output', required=True, help='output TSV file')
    return parser.parse_args()


HEADER = (
    "sequence_id\tproductive\ttop_v_allele\ttop_v_identity\tmutated_status\t"
    "top_v_bit_score\ttop_v_evalue\ttop_j_allele\ttop_v_alleles\t"
    "all_v_alleles\tall_v_identities\tall_mutated_status\tall_bit_scores\tall_evalues\tbtop_mutations\n"
)


def write_query(out, seq_id, productive, top_v_alleles_summary, j_allele,
                v_ids, v_pcts, v_bits, v_evals, v_btops):
    if not v_ids:
        out.write(
            f"{seq_id}\t{productive}\tN/A\tN/A\tN/A\tN/A\tN/A\t{j_allele}\t"
            "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n"
        )
        return

    # Select best hit: highest bit score, then lowest e-value
    if len(v_ids) == 1:
        top_idx = 0
    elif "," in top_v_alleles_summary:
        top_idx = sorted(
            range(len(v_ids)),
            key=lambda i: (-float(v_bits[i]), float(v_evals[i]))
        )[0]
    else:
        try:
            top_idx = v_ids.index(top_v_alleles_summary)
        except ValueError:
            top_idx = 0

    top_v = v_ids[top_idx]
    top_pct = v_pcts[top_idx]
    top_bit = v_bits[top_idx]
    top_eval = v_evals[top_idx]
    mutated = "MUTATED" if float(top_pct) < 98 else "UNMUTATED"
    all_mutated = ",".join("MUTATED" if float(x) < 98 else "UNMUTATED" for x in v_pcts)

    out.write(
        f"{seq_id}\t{productive}\t{top_v}\t{top_pct}\t{mutated}\t"
        f"{top_bit}\t{top_eval}\t{j_allele}\t{top_v_alleles_summary}\t"
        f"{','.join(v_ids)}\t{','.join(v_pcts)}\t{all_mutated}\t"
        f"{','.join(v_bits)}\t{','.join(v_evals)}\t{','.join(v_btops)}\n"
    )


def main(args):
    # Per-query state
    seq_id = None
    productive = "N/A"
    top_v_alleles_summary = "N/A"
    j_allele = "N/A"
    v_allele_pos = None
    productive_pos = None
    j_allele_pos = None
    vdj_summary_next = False
    hit_positions = {}
    v_ids, v_pcts, v_bits, v_evals, v_btops = [], [], [], [], []

    with open(args.input, 'r') as handle, open(args.output, 'w') as out:
        out.write(HEADER)

        for line_num, raw_line in enumerate(handle, start=1):
            if line_num == 1:
                continue  # skip file header "# IGBLASTN x.x.x"

            line = raw_line.rstrip('\n')

            # Start of a new query (or the very first query)
            if line.startswith("# IGBLASTN") or line.startswith("# Query:"):
                # Flush previous query
                if seq_id is not None:
                    write_query(out, seq_id, productive, top_v_alleles_summary, j_allele,
                                v_ids, v_pcts, v_bits, v_evals, v_btops)

                if line.startswith("# Query:"):
                    seq_id = line[9:].strip()
                    productive = "N/A"
                    top_v_alleles_summary = "N/A"
                    j_allele = "N/A"
                    v_allele_pos = None
                    productive_pos = None
                    j_allele_pos = None
                    vdj_summary_next = False
                    hit_positions = {}
                    v_ids, v_pcts, v_bits, v_evals, v_btops = [], [], [], [], []
                else:
                    # "# IGBLASTN" between queries resets seq_id so next "# Query:" initialises
                    seq_id = None
                continue

            if seq_id is None:
                continue

            # VDJ rearrangement summary header — column positions vary by chain type
            if "rearrangement summary" in line:
                info_line = (
                    line.strip()
                    .rstrip(".  Multiple equivalent top matches, if present, are separated by a comma.")
                    .split(" (")
                )
                if len(info_line) > 1:
                    cols = info_line[1].rstrip(")").split(", ")
                    for pos, col in enumerate(cols):
                        if col.startswith("Top V gene match"):
                            v_allele_pos = pos
                        if col.startswith("Productive"):
                            productive_pos = pos
                        if col.startswith("Top J gene match"):
                            j_allele_pos = pos
                vdj_summary_next = True
                continue

            # Line immediately after the VDJ summary header is the data row
            if vdj_summary_next:
                vdj_summary_next = False
                fields = line.split('\t')
                if v_allele_pos is not None and v_allele_pos < len(fields):
                    top_v_alleles_summary = fields[v_allele_pos]
                if productive_pos is not None and productive_pos < len(fields):
                    productive = fields[productive_pos]
                if j_allele_pos is not None and j_allele_pos < len(fields):
                    j_allele = fields[j_allele_pos]
                continue

            # Hit table column header
            if line.startswith("# Fields:"):
                cols = line.lstrip("# Fields: ").split(", ")
                # First column in each data row is segment type, not in the Fields list
                for pos, col in enumerate(cols):
                    hit_positions[col] = pos + 1
                continue

            # V-gene hit rows
            if line.startswith("V\t") and hit_positions:
                fields = line.split('\t')
                v_ids.append(fields[hit_positions["subject id"]])
                v_pcts.append(fields[hit_positions["% identity"]])
                v_bits.append(fields[hit_positions["bit score"]])
                v_evals.append(fields[hit_positions["evalue"]])
                v_btops.append(fields[hit_positions["BTOP"]])
                continue

            # End of file sentinel
            if line.startswith("Total queries"):
                if seq_id is not None:
                    write_query(out, seq_id, productive, top_v_alleles_summary, j_allele,
                                v_ids, v_pcts, v_bits, v_evals, v_btops)
                break


if __name__ == "__main__":
    args = parse_args()
    main(args)
