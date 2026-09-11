#!/usr/bin/env python3
"""Join pvactools' own aggregated neoantigen report onto a native-build MAF, producing a small,
MAF-shaped table (the MAF's own first N columns plus a handful of pvactools annotation columns) --
one row per pvactools-predicted mutation, not a left join against every MAF row.

The native MAF here is this module's own existing inputs.vcf2maf_maf (the same one
_pvactools_extract_coding_regions already reads) -- both it and this aggregated report trace back
to the exact same VCF records, so a genomic-coordinate join between them is unambiguous (no
cross-tool/cross-build annotation-consistency assumption needed) -- see this module's CHANGELOG for
the full reasoning behind this design, including the approaches that were tried and rejected first.

pvactools' own aggregated report has no separate Chromosome/Start/Reference/Variant columns --
confirmed against real production files that its "ID" column encodes them as
"{Chromosome}-{Start}-{Stop}-{Reference}-{Variant}" (e.g. "19-11134251-11134252-G-A"). Indels use
VCF-style REF/ALT (the anchor base included, e.g. Reference="GTGGCCACGATGCGC" Variant="G"), never a
bare "-" placeholder -- also confirmed against real production indel rows -- so a regex anchored on
Start/Stop being pure digits and Reference/Variant being pure nucleotide sequence is safe.
"""
import argparse
import csv
import json
import re

parser = argparse.ArgumentParser()
parser.add_argument("--aggregated-report", required=True, help="pvactools' own Combined.all_epitopes.aggregated.tsv")
parser.add_argument("--native-maf", required=True, help="This module's own inputs.vcf2maf_maf (native-build MAF)")
parser.add_argument("--relaxed-filtered", required=True, help="This pair's additional_filters/relaxed/*.filtered.tsv")
parser.add_argument("--very-relaxed-filtered", required=True, help="This pair's additional_filters/very_relaxed/*.filtered.tsv")
parser.add_argument("--minimal-maf-columns", type=int, required=True, help="Keep only the first N columns of the native MAF (options.minimal_maf_columns)")
parser.add_argument("--neoantigen-columns", nargs="+", required=True, help="Aggregated-report column names to carry through (options.neoantigen_maf_columns)")
parser.add_argument("--pass-filters-json", default="{}", help="JSON dict: {column_name: {binding_threshold, presentation_percentile_threshold, min_dna_vaf, excluded_tiers}} (options.neoantigen_pass_filters)")
parser.add_argument("--output-maf", required=True)
parser.add_argument("--output-audit", required=True)
args = parser.parse_args()
pass_filters = json.loads(args.pass_filters_json)

# Chromosome is "everything before the first digit-only run" -- greedy but safe, since standard
# human chromosome names (1-22, X, Y, MT) never contain a "-digits-digits-" pattern themselves.
ID_PATTERN = re.compile(r"^(.+)-(\d+)-(\d+)-([ACGTN]+)-([ACGTN]+)$")

# Tolerant fallback matching window (bp) -- mirrors src/maf_coding_positions.py's own +/-1bp
# tolerance for MAF/VCF indel-representation differences, widened slightly since here we're
# matching two independently-generated MAF-shaped representations of the same VCF record rather
# than a MAF against its own source VCF directly.
TOLERANCE_BP = 2


def parse_locus_from_id(id_value):
    match = ID_PATTERN.match(id_value)
    if not match:
        return None
    chrom, start, _stop, ref, var = match.groups()
    return chrom, start, ref, var


def read_native_maf(path, n_columns):
    with open(path) as f:
        # MAFs conventionally carry a variable number of leading '#'-prefixed comment lines before
        # the real header row -- same convention already handled in maf_coding_positions.py.
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        minimal_header = header[:n_columns]
        reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
        exact_index = {}
        ref_var_index = {}
        by_chrom = {}
        for row in reader:
            chrom, start, ref, var = row["Chromosome"], row["Start_Position"], row["Reference_Allele"], row["Tumor_Seq_Allele2"]
            minimal_row = {col: row[col] for col in minimal_header}
            start_int = int(start)
            end_int = start_int + max(len(ref), 1) - 1
            exact_index[(chrom, start, ref, var)] = minimal_row
            ref_var_index.setdefault((chrom, ref, var), []).append((start_int, minimal_row))
            by_chrom.setdefault(chrom, []).append((start_int, end_int, minimal_row))
    return minimal_header, exact_index, ref_var_index, by_chrom


def find_maf_match(chrom, start, ref, var, exact_index, ref_var_index, by_chrom):
    key = (chrom, start, ref, var)
    if key in exact_index:
        return exact_index[key], "exact"

    # Tier 2: same Reference/Variant strings, position drifted slightly (e.g. one representation
    # trimmed/shifted an anchor base but otherwise agrees on the actual allele sequences).
    start_int = int(start)
    best = None
    for candidate_start, row in ref_var_index.get((chrom, ref, var), []):
        offset = abs(candidate_start - start_int)
        if offset <= TOLERANCE_BP and (best is None or offset < best[0]):
            best = (offset, row)
    if best is not None:
        return best[1], f"tolerant, same alleles (+/-{best[0]}bp)"

    # Tier 3: position-range overlap only, ignoring exact Reference/Variant string equality --
    # real MAF/VCF indel representations can shift which anchor base is included, changing the
    # allele strings themselves, not just the position (confirmed necessary via a real synthetic
    # test during this feature's own development, not a hypothetical). Mirrors
    # src/maf_coding_positions.py's own position-only matching philosophy for indels, and accepts
    # the same known risk that module already documents: at a genuinely multiallelic/decomposed
    # site this could match the wrong nearby candidate -- logged via the audit trail's own
    # match-type label, not silently indistinguishable from an exact match.
    end_int = start_int + max(len(ref), 1) - 1
    best = None
    for candidate_start, candidate_end, row in by_chrom.get(chrom, []):
        if candidate_start <= end_int + TOLERANCE_BP and start_int <= candidate_end + TOLERANCE_BP:
            offset = abs(candidate_start - start_int)
            if best is None or offset < best[0]:
                best = (offset, row)
    if best is not None:
        return best[1], f"tolerant, position-overlap only (+/-{best[0]}bp)"

    return None, None


def _passes_numeric_criterion(raw_value, threshold, comparison):
    # "NA"/missing auto-passes a criterion -- matches pvactools' own established convention
    # (confirmed directly from lib/filter.py's Filter.execute() earlier this session: a column
    # value of 'NA' never causes exclusion there either), kept consistent here rather than
    # inventing a stricter rule this module doesn't use anywhere else.
    if threshold is None:
        return True
    if raw_value is None or raw_value == "" or raw_value == "NA":
        return True
    value = float(raw_value)
    return comparison(value, threshold)


def passes_filter(aggregated_row, spec):
    # aggregated_row is the RAW pvactools aggregated-report row (has 'Tier'/'IC50 MT'/
    # 'Pres %ile MT'/'DNA VAF' keys) -- evaluated against the fixed criteria this module supports,
    # each individually optional (omit a key in the preset's own spec to skip that criterion
    # entirely, not just to make it always-pass). A criterion is a real numeric comparison, not
    # "NA-tolerant" in the sense of ignoring bad data -- only genuinely missing ('NA') values
    # auto-pass, matching pvactools' own convention.
    if not _passes_numeric_criterion(aggregated_row.get("IC50 MT"), spec.get("binding_threshold"), lambda v, t: v <= t):
        return False
    if not _passes_numeric_criterion(aggregated_row.get("Pres %ile MT"), spec.get("presentation_percentile_threshold"), lambda v, t: v <= t):
        return False
    if not _passes_numeric_criterion(aggregated_row.get("DNA VAF"), spec.get("min_dna_vaf"), lambda v, t: v >= t):
        return False
    excluded_tiers = spec.get("excluded_tiers") or []
    if excluded_tiers and aggregated_row.get("Tier") in excluded_tiers:
        return False
    return True


def load_filter_membership(path):
    # additional_filters/{preset}/*.filtered.tsv is derived from Combined.all_epitopes.tsv (the
    # uncollapsed report), which does carry real Chromosome/Start/Reference/Variant columns --
    # confirmed against real production files, unlike the aggregated report's own encoded ID.
    membership = set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            return membership
        for row in reader:
            if {"Chromosome", "Start", "Reference", "Variant"} <= row.keys():
                membership.add((row["Chromosome"], row["Start"], row["Reference"], row["Variant"]))
    return membership


minimal_header, exact_index, ref_var_index, by_chrom = read_native_maf(args.native_maf, args.minimal_maf_columns)
relaxed_pass = load_filter_membership(args.relaxed_filtered)
very_relaxed_pass = load_filter_membership(args.very_relaxed_filtered)

audit = []  # (aggregated_report_id, status)
output_rows = []
n_total = 0

with open(args.aggregated_report) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        n_total += 1
        id_value = row["ID"]
        locus = parse_locus_from_id(id_value)
        if locus is None:
            audit.append((id_value, "skipped: could not parse genomic locus from ID"))
            continue
        chrom, start, ref, var = locus
        maf_row, match_type = find_maf_match(chrom, start, ref, var, exact_index, ref_var_index, by_chrom)
        if maf_row is None:
            audit.append((id_value, "skipped: no matching row in native-build MAF"))
            continue
        out_row = dict(maf_row)
        for col in args.neoantigen_columns:
            out_row[col] = row.get(col, "NA")
        out_row["Pvacseq_Pass_Relaxed"] = str((chrom, start, ref, var) in relaxed_pass)
        out_row["Pvacseq_Pass_VeryRelaxed"] = str((chrom, start, ref, var) in very_relaxed_pass)
        for column_name, spec in pass_filters.items():
            out_row[column_name] = str(passes_filter(row, spec))
        audit.append((id_value, f"matched ({match_type})"))
        output_rows.append(out_row)

output_fieldnames = (
    minimal_header + list(args.neoantigen_columns) +
    ["Pvacseq_Pass_Relaxed", "Pvacseq_Pass_VeryRelaxed"] + list(pass_filters.keys())
)
with open(args.output_maf, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=output_fieldnames, delimiter="\t", restval="NA")
    writer.writeheader()
    writer.writerows(output_rows)

with open(args.output_audit, "w") as f:
    f.write("aggregated_report_id\tstatus\n")
    for id_value, status in audit:
        f.write(f"{id_value}\t{status}\n")

print(f"Matched {len(output_rows)} of {n_total} pvactools-predicted mutations to the native-build MAF.")
