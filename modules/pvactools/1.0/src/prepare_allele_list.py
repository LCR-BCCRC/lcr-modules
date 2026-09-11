#!/usr/bin/env python3
"""Build pvacseq's comma-separated HLA allele-list argument entirely from mhc_hammer's own
HLA-HD-based typing: Class I (_mhc_hammer_output_hla_final_result, HLA-HD's own consolidated
per-patient result) and Class II (_mhc_hammer_output_hla2_alleles).

Real formats confirmed against a real installed pvactools=7.1.3 (`pvactools valid_alleles`):
  Class I:  "HLA-A*02:01"                (gene short name, with HLA- prefix)
  Class II: "DRB1*01:01"                 (single-chain genes, no HLA- prefix)
            "DQA1*01:01-DQB1*02:01"      (DQ heterodimer, hyphenated pair)
            "DPA1*01:03-DPB1*01:01"      (DP heterodimer, hyphenated pair)

Class I input is HLA-HD's own `{patient_id}_final.result.txt` (gene name column plus two allele
columns, alleles carrying an "HLA-" prefix). Confirmed against real production files (see this
module's CHANGELOG): a "Not typed" gene shows that literal placeholder in both columns, while a
gene with only one confidently-called allele (homozygous, or a second allele HLA-HD couldn't
resolve) shows a literal "-" in the second column instead -- distinct from "Not typed" and handled
separately below. Parsing is deliberately defensive beyond those two known cases too (skips any
value with no "*" at all, rather than assuming an exhaustive list of placeholder strings) so a
real file using some other not-yet-seen HLA-HD convention still degrades gracefully instead of
crashing.
"""
import argparse
import csv
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("--hla-final-result", required=True, help="mhc_hammer's own {patient_id}_final.result.txt (Class I)")
parser.add_argument("--hla2-alleles", default=None, help="mhc_hammer's own hla2_alleles.csv (optional, Class II)")
parser.add_argument("--output-allele-list", required=True)
parser.add_argument("--output-audit", required=True)
args = parser.parse_args()

audit = []  # (source, gene, raw, resolved, note)

CLASS_I_GENES = {"A", "B", "C"}


def trim_resolution(allele):
    # GENE*NN:NN[:NN[:NN]] -> GENE*NN:NN.
    star = allele.index("*")
    gene, fields = allele[:star], allele[star + 1:].split(":")
    return f"{gene}*{fields[0]}:{fields[1]}" if len(fields) >= 2 else allele


# ---- Class I (mhc_hammer's own HLA-HD final result) ----
class_i = set()
with open(args.hla_final_result) as f:
    for row in csv.reader(f, delimiter="\t"):
        if not row or not row[0]:
            continue
        gene, alleles = row[0].strip(), row[1:]
        if gene not in CLASS_I_GENES:
            continue
        for allele in alleles:
            allele = allele.strip()
            if not allele or allele.lower() == "not typed":
                audit.append(("hla_final_result", gene, allele, "", "skipped: not typed"))
                continue
            if allele == "-":
                # HLA-HD's own convention for "no distinct second allele" (homozygous, or only
                # one allele group reported) -- confirmed via a real production sample
                # (modules/pvactools/CHANGELOG.md) whose *_final.result.txt had a literal "-" in
                # a gene's second allele column. Matches how modules/mhc_hammer/1.0's own Class II
                # parser (src/parse_hlahd_output.R) already handles the identical HLA-HD output
                # convention (duplicates allele1 rather than treating it as a no-call). The gene's
                # one real call is already captured from its other column, so there's nothing to
                # add here -- class_i is a set, so explicitly duplicating allele1 would be a no-op
                # anyway.
                audit.append(("hla_final_result", gene, allele, "", "skipped: homozygous/single call (HLA-HD '-')"))
                continue
            if "*" not in allele:
                # Defensive catch-all: an HLA-HD placeholder we haven't seen before. Skip rather
                # than crash trim_resolution() -- matches this module's own stated design
                # philosophy (see this file's header docstring) of degrading gracefully instead of
                # failing the whole pair over one unrecognized value.
                audit.append(("hla_final_result", gene, allele, "", f"skipped: unrecognized allele format (raw: {allele!r})"))
                continue
            # HLA-HD's final result conventionally already carries the HLA- prefix -- strip it
            # first so it's never doubled, then add it back consistently.
            bare = allele[4:] if allele.upper().startswith("HLA-") else allele
            resolved = "HLA-" + trim_resolution(bare)
            class_i.add(resolved)
            audit.append(("hla_final_result", gene, allele, resolved, ""))

# ---- Class II (mhc_hammer HLA-HD), optional ----
class_ii_single = set()   # DRB1/DRB3/DRB4/DRB5
dq_a, dq_b = set(), set()
dp_a, dp_b = set(), set()
# HLA-DM/-DO are non-classical peptide-loading chaperones, never presented to T cells -- not
# relevant to neoantigen prediction, and pvacseq has no alleles for them at all. DRA is DRB1's
# invariant alpha chain (unlike DQ/DP, DR's alpha chain is essentially non-polymorphic and is
# never submitted as its own standalone allele) -- confirmed absent from every one of pvactools'
# own Class II allele reference files (MHCnuggets.txt, netmhciipan.tsv, nn_align.tsv,
# smm_align.tsv, MixMHC2pred.tsv). Before this fix, "DRA*01:01"/"DRA*01:02" were passed to
# pvacseq as if they were real single-chain alleles like DRB1 -- silently dropped per-run with a
# "not valid for Method ..." warning (pvacseq's per-allele/method check just skips unrecognized
# alleles rather than failing the run), so this was wasted work/log noise rather than a
# correctness bug, but excluding it here matches how DMA/DMB/DOA/DOB are already handled.
IGNORED_GENES = {"DMA", "DMB", "DOA", "DOB", "DRA"}

if args.hla2_alleles:
    with open(args.hla2_alleles) as f:
        for gene, allele1, allele2 in csv.reader(f):
            if gene in IGNORED_GENES:
                audit.append(("hla2_alleles", gene, f"{allele1},{allele2}", "",
                              "skipped: non-classical, not presented to T cells"))
                continue
            for allele in (allele1, allele2):
                if allele == "not typed":
                    audit.append(("hla2_alleles", gene, allele, "", "skipped: not typed"))
                    continue
                # HLA-HD's raw candidate string may or may not already carry the gene prefix --
                # not independently verified (see this module's CHANGELOG), so add it defensively.
                raw_with_gene = allele if allele.startswith(gene) else f"{gene}*{allele.split('*')[-1]}"
                resolved = trim_resolution(raw_with_gene)
                audit.append(("hla2_alleles", gene, allele, resolved, ""))
                if gene in ("DQA1",):
                    dq_a.add(resolved)
                elif gene in ("DQB1",):
                    dq_b.add(resolved)
                elif gene in ("DPA1",):
                    dp_a.add(resolved)
                elif gene in ("DPB1",):
                    dp_b.add(resolved)
                else:
                    class_ii_single.add(resolved)

class_ii = set(class_ii_single)
# Can't phase which chain pairs with which from unphased WES/WGS -- submit every combination
# rather than guessing (see this module's plan/CHANGELOG for the reasoning).
class_ii |= {f"{a}-{b}" for a, b in itertools.product(dq_a, dq_b)}
class_ii |= {f"{a}-{b}" for a, b in itertools.product(dp_a, dp_b)}

all_alleles = sorted(class_i) + sorted(class_ii)
if not all_alleles:
    raise AssertionError(
        "No usable HLA alleles from mhc_hammer's own Class I or Class II typing -- cannot run "
        "pvacseq for this pair."
    )

with open(args.output_allele_list, "w") as f:
    f.write(",".join(all_alleles))

with open(args.output_audit, "w") as f:
    f.write("source\tgene\traw\tresolved\tnote\n")
    for row in audit:
        f.write("\t".join(row) + "\n")

print(f"Class I: {len(class_i)} allele(s). Class II: {len(class_ii)} allele(s)/pair(s).")
