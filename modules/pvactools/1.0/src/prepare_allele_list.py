#!/usr/bin/env python3
"""Build pvacseq's comma-separated HLA allele-list argument from LILAC (Class I) + mhc_hammer's
HLA-HD-based Class II typing.

Real formats confirmed against a real installed pvactools=7.1.3 (`pvactools valid_alleles`):
  Class I:  "HLA-A*02:01"                (gene short name, with HLA- prefix)
  Class II: "DRB1*01:01"                 (single-chain genes, no HLA- prefix)
            "DQA1*01:01-DQB1*02:01"      (DQ heterodimer, hyphenated pair)
            "DPA1*01:03-DPB1*01:01"      (DP heterodimer, hyphenated pair)
"""
import argparse
import csv
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("--lilac-tsv", required=True)
parser.add_argument("--hla2-alleles", default=None, help="mhc_hammer's own hla2_alleles.csv (optional)")
parser.add_argument("--output-allele-list", required=True)
parser.add_argument("--output-audit", required=True)
args = parser.parse_args()

audit = []  # (source, gene, raw, resolved, note)


def trim_resolution(allele):
    # GENE*NN:NN[:NN[:NN]] -> GENE*NN:NN. LILAC already reports at 2-field resolution; HLA-HD's
    # own raw candidate strings are not independently verified to (see this module's CHANGELOG).
    star = allele.index("*")
    gene, fields = allele[:star], allele[star + 1:].split(":")
    return f"{gene}*{fields[0]}:{fields[1]}" if len(fields) >= 2 else allele


# ---- Class I (LILAC) ----
class_i = set()
with open(args.lilac_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        raw = row["Allele"]
        resolved = "HLA-" + trim_resolution(raw)
        class_i.add(resolved)
        audit.append(("lilac", raw.split("*")[0], raw, resolved, ""))

# ---- Class II (mhc_hammer HLA-HD), optional ----
class_ii_single = set()   # DRB1/DRB3/DRB4/DRB5
dq_a, dq_b = set(), set()
dp_a, dp_b = set(), set()
# HLA-DM/-DO are non-classical peptide-loading chaperones, never presented to T cells -- not
# relevant to neoantigen prediction, and pvacseq has no alleles for them at all.
IGNORED_GENES = {"DMA", "DMB", "DOA", "DOB"}

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
        "No usable HLA alleles from either LILAC or mhc_hammer -- cannot run pvacseq for this pair."
    )

with open(args.output_allele_list, "w") as f:
    f.write(",".join(all_alleles))

with open(args.output_audit, "w") as f:
    f.write("source\tgene\traw\tresolved\tnote\n")
    for row in audit:
        f.write("\t".join(row) + "\n")

print(f"Class I: {len(class_i)} allele(s). Class II: {len(class_ii)} allele(s)/pair(s).")
