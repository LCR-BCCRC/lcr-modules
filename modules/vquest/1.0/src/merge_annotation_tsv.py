#!/usr/bin/env python3
"""
Left-join an annotation TSV onto a base TSV by a shared key column.
Rows in the base table that have no match in the annotation table get 'N/A'
for every annotation column.

When --base_key and --annot_key differ, the base key value is passed through
_extract_mixcr_clone_id before the annotation lookup. This handles the case
where the base uses V-QUEST sequence_id (which for MiXCR inputs encodes the
cloneId in compound format "cloneId_N_readFraction_F_readCount_C") while the
annotation uses the bare cloneId column. For non-MiXCR sequence_ids the
function is a no-op.

When --sample_id is supplied, sample_id is injected as the first column (if
absent) and all columns are reordered into a harmonised schema that is
consistent across both igseqr and MiXCR inputs.  Two additional
normalisation steps are applied when --sample_id is given:

  * igseqr inputs: the V-QUEST 'sequence_id' column (which carries the
    transcript_id value) is renamed 'transcript_id' for consistency with
    igblast outputs.
  * MiXCR inputs: a 'cloneId' column is extracted from the compound
    'sequence_id' header written by mixcr_to_fasta.py, so the primary ID
    column matches the bare integer used in igblast outputs.
"""

import argparse
import csv
import re
import sys


def _extract_mixcr_clone_id(sequence_id):
    """Parse cloneId from a MiXCR compound sequence_id.

    mixcr_to_fasta.py writes FASTA headers as:
        cloneId_{N}_readFraction_{F}_readCount_{C}
    Returns the cloneId string, or the original value if the format does not match.
    """
    m = re.match(r"cloneId_(\S+?)_readFraction_", sequence_id)
    return m.group(1) if m else sequence_id


def _ordered_columns(all_cols):
    """Return de-duplicated columns in harmonised output order.

    Columns not explicitly listed fall through to the end in their original
    relative order.  Handles the V-QUEST-specific L-gene columns and the
    full LVDJC alignment blocks.
    """
    col_set = set(all_cols)

    def present(*names):
        return [c for c in names if c in col_set]

    # Detection: igseqr adds kallisto/salmon quant cols; MiXCR adds readCount.
    is_igseqr = any(c in col_set for c in ("length", "eff_length", "est_counts", "tpm"))
    is_mixcr  = "readCount" in col_set

    order = ["sample_id"]

    # --- ID column ---
    if is_igseqr:
        # sequence_id has been renamed to transcript_id before this call
        order += present("transcript_id")
    elif is_mixcr:
        # cloneId was extracted from compound sequence_id before this call
        order += present("cloneId")
    else:
        order += present("sequence_id")

    # --- Quantification ---
    if is_igseqr:
        order += present("length", "eff_length", "est_counts", "tpm")
    elif is_mixcr:
        order += present("readCount", "readFraction")

    # --- Annotation quality flags ---
    order += present("locus", "stop_codon", "vj_in_frame", "productive", "rev_comp",
                     "complete_vdj", "complete_lvdjc")

    # --- Gene calls (L first for V-QUEST), then MiXCR best calls ---
    order += present("l_call", "v_call", "d_call", "j_call", "c_call")
    if is_mixcr:
        order += present("bestVHit", "bestDHit", "bestJHit", "bestCHit",
                         "allVHitsWithScore", "allDHitsWithScore",
                         "allJHitsWithScore", "allCHitsWithScore")

    # --- Glycosylation ---
    order += present("num_glycosylation_sites", "glycosylation_imgt_positions",
                     "glycosylation_motifs",
                     "num_acquired_glycosylation_sites", "acquired_glycosylation_imgt_positions",
                     "germline_aa_sequence",
                     "num_germline_glycosylation_sites", "germline_glycosylation_imgt_positions",
                     "germline_glycosylation_motifs")

    # --- Alignment quality: score, support, identity, cigar + positions — L(V)DJC ---
    # L gene (no support field in V-QUEST)
    order += present("l_score", "l_identity", "l_cigar",
                     "l_sequence_start", "l_sequence_end",
                     "l_germline_start", "l_germline_end")
    for gene in ("v", "d", "j", "c"):
        order += present(
            f"{gene}_score", f"{gene}_support", f"{gene}_identity", f"{gene}_cigar",
            f"{gene}_sequence_start", f"{gene}_sequence_end",
            f"{gene}_germline_start", f"{gene}_germline_end",
            f"{gene}_alignment_start", f"{gene}_alignment_end",
        )

    # --- Junction / NP regions (V-QUEST adds np*_aa and extra length fields) ---
    order += present("np1", "np1_aa", "np1_length",
                     "np2", "np2_aa", "np2_length",
                     "n1_length", "n2_length",
                     "p3v_length", "p5d_length", "p3d_length", "p5j_length",
                     "junction", "junction_length", "junction_aa", "junction_aa_length",
                     "junction_decryption")

    # --- Region boundary positions ---
    for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4", "cdr3"):
        order += present(f"{region}_start", f"{region}_end")

    # --- AIRR metadata fields ---
    order += present("consensus_count", "duplicate_count", "cell_id", "clone_id",
                     "rearrangement_id", "repertoire_id", "rearrangement_set_id",
                     "sequence_analysis_category", "d_number",
                     "5prime_trimmed_n_nb", "3prime_trimmed_n_nb",
                     "insertions", "deletions")

    # --- Sequences: full assembled sequence first (LVDJC and VDJ alignments) ---
    order += present("sequence", "sequence_aa", "aa_sequence",
                     "LVDJC_sequence_alignment", "LVDJC_sequence_alignment_aa",
                     "LVDJC_germline_alignment", "LVDJC_germline_alignment_aa",
                     "LVJC_sequence_alignment", "LVJC_sequence_alignment_aa",
                     "LVJC_germline_alignment", "LVJC_germline_alignment_aa",
                     "sequence_alignment", "sequence_alignment_aa",
                     "germline_alignment", "germline_alignment_aa")

    # --- Sequences: gene components L(V)DJC ---
    for gene in ("l", "v", "d", "j", "c"):
        order += present(
            f"{gene}_sequence_alignment", f"{gene}_sequence_alignment_aa",
            f"{gene}_germline_alignment", f"{gene}_germline_alignment_aa",
        )

    # --- Region sequences ---
    for region in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4", "cdr3"):
        order += present(region, f"{region}_aa")

    # --- Remaining MiXCR-specific columns (alignments, quality, region seqs) ---
    if is_mixcr:
        order += present(
            "allVAlignments", "allDAlignments", "allJAlignments", "allCAlignments",
            "targetSequences", "targetQualities",
            "nSeqFR1", "minQualFR1", "nSeqCDR1", "minQualCDR1",
            "nSeqFR2", "minQualFR2", "nSeqCDR2", "minQualCDR2",
            "nSeqFR3", "minQualFR3", "nSeqCDR3", "minQualCDR3",
            "nSeqFR4", "minQualFR4",
            "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2",
            "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4", "refPoints",
        )

    # --- Catch-all: anything not yet placed, preserving relative input order ---
    placed = set(order)
    for c in all_cols:
        if c not in placed:
            order.append(c)
            placed.add(c)

    return order


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base",       required=True, help="base TSV (all rows kept)")
    parser.add_argument("--annotation", required=True, help="annotation TSV to join in")
    parser.add_argument("--base_key",   default="sequence_id", help="join key column in base")
    parser.add_argument("--annot_key",  default="sequence_id", help="join key column in annotation")
    parser.add_argument("--output",     required=True, help="output TSV path")
    parser.add_argument("--sample_id",  default=None,
                        help="Sample ID wildcard value; triggers column reordering when supplied")
    args = parser.parse_args()

    # When key columns differ, transform the base key value at lookup time so
    # a compound V-QUEST sequence_id resolves to the annotation's bare cloneId.
    transform_base_key = args.base_key != args.annot_key

    # Read annotation fully, keeping all non-key columns.
    annot, all_annot_cols = {}, []
    with open(args.annotation) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if args.annot_key not in reader.fieldnames:
            sys.exit(
                f"ERROR: --annot_key '{args.annot_key}' not found in {args.annotation}.\n"
                f"       Available columns: {', '.join(reader.fieldnames)}\n"
                f"       Check options.source_id_column in the vquest module config."
            )
        all_annot_cols = [c for c in reader.fieldnames if c != args.annot_key]
        for row in reader:
            annot[row[args.annot_key]] = {c: row[c] for c in all_annot_cols}

    with open(args.base) as fh_in, open(args.output, "w", newline="") as fh_out:
        reader = csv.DictReader(fh_in, delimiter="\t")

        # Only add annotation columns not already present in the base, preventing
        # duplicate column names (e.g. 'sequence' exists in both igseqr source_tsv
        # and the V-QUEST AIRR output).
        base_col_set = set(reader.fieldnames)
        annot_cols = [c for c in all_annot_cols if c not in base_col_set]

        raw_header = list(reader.fieldnames) + annot_cols

        # Flags computed once before the row loop for use inside it.
        rename_seq_id = False
        add_clone_id  = False

        if args.sample_id is not None:
            if "sample_id" not in raw_header:
                raw_header = ["sample_id"] + raw_header

            is_igseqr = any(c in raw_header for c in ("length", "eff_length", "est_counts", "tpm"))
            is_mixcr  = "readCount" in raw_header

            # igseqr: rename 'sequence_id' → 'transcript_id' so outputs from
            # both igblast and V-QUEST use the same primary ID column name.
            if (is_igseqr
                    and "sequence_id" in raw_header
                    and "transcript_id" not in raw_header):
                raw_header = ["transcript_id" if c == "sequence_id" else c
                              for c in raw_header]
                rename_seq_id = True

            # MiXCR: extract the bare cloneId from the compound sequence_id so
            # the primary ID column matches the integer used in igblast outputs.
            if (is_mixcr
                    and "sequence_id" in raw_header
                    and "cloneId" not in raw_header):
                raw_header = ["cloneId"] + raw_header
                add_clone_id = True

            header = _ordered_columns(raw_header)
        else:
            header = raw_header

        writer = csv.DictWriter(fh_out, fieldnames=header, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            lookup = (_extract_mixcr_clone_id(row[args.base_key])
                      if transform_base_key else row[args.base_key])
            matched = annot.get(lookup)
            row.update(
                {c: matched[c] for c in annot_cols}
                if matched is not None
                else {c: "N/A" for c in annot_cols}
            )
            if args.sample_id is not None:
                row.setdefault("sample_id", args.sample_id)
                if rename_seq_id:
                    row["transcript_id"] = row.pop("sequence_id", "")
                if add_clone_id:
                    row["cloneId"] = _extract_mixcr_clone_id(row.get("sequence_id", ""))
            writer.writerow(row)


if __name__ == "__main__":
    main()
