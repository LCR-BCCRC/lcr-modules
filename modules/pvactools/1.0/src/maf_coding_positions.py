#!/usr/bin/env python3
"""Extract bcftools regions (+ audit TSV) for protein-altering variants from a vcf2maf MAF.

Reuses vcf2maf's own already-computed Variant_Classification (from its earlier VEP run) so
pvactools never re-annotates/re-scores variants already known to be silent/non-coding -- the
large majority for a WGS sample. Position-only matching (Start_Position/End_Position, +/-1bp
padding), not exact ref/alt, to tolerate MAF/VCF indel-representation differences.
"""
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("--maf", required=True, help="vcf2maf's own final MAF (_vcf2maf_output_maf)")
parser.add_argument("--classifications", nargs="+", required=True,
                     help="Variant_Classification values to keep (options.coding_variant_classifications)")
parser.add_argument("--output-regions", required=True, help="bcftools -R regions file (CHROM/FROM/TO)")
parser.add_argument("--output-audit", required=True, help="Per-classification kept/dropped counts TSV")
args = parser.parse_args()

keep = set(args.classifications)
counts = {}
regions = []

with open(args.maf) as f:
    # MAFs conventionally carry a variable number of leading '#'-prefixed comment lines before
    # the real header row -- skip them the same way any standard MAF reader would.
    header = None
    for line in f:
        if line.startswith("#"):
            continue
        header = line.rstrip("\n").split("\t")
        break
    reader = csv.DictReader(f, fieldnames=header, delimiter="\t")
    for row in reader:
        classification = row["Variant_Classification"]
        counts[classification] = counts.get(classification, [0, 0])
        counts[classification][0] += 1
        if classification not in keep:
            continue
        counts[classification][1] += 1
        chrom = row["Chromosome"]
        start = int(row["Start_Position"]) - 1
        end = int(row["End_Position"]) + 1
        regions.append((chrom, max(start, 1), end))

with open(args.output_regions, "w") as f:
    for chrom, start, end in regions:
        f.write(f"{chrom}\t{start}\t{end}\n")

with open(args.output_audit, "w") as f:
    f.write("Variant_Classification\ttotal\tkept\n")
    for classification, (total, kept) in sorted(counts.items()):
        f.write(f"{classification}\t{total}\t{kept}\n")

print(f"Kept {len(regions)} of {sum(t for t, _ in counts.values())} variants "
      f"({', '.join(sorted(keep))})")
