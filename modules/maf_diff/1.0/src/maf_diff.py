import sys
import pandas as pd

log = open(snakemake.log[0], "w")
sys.stderr = log

caller1_name = snakemake.params.caller1_name
caller2_name = snakemake.params.caller2_name
key_cols = snakemake.params.key_cols
driver_genes_path = snakemake.params.driver_genes
driver_gene_col = snakemake.params.driver_gene_col
driver_col = snakemake.params.driver_col
driver_col_value = snakemake.params.driver_col_value
maf_gene_col = snakemake.params.maf_gene_col
maf_vc_col = snakemake.params.maf_vc_col
coding_variant_classes = set(snakemake.params.coding_variant_classes)

maf1 = pd.read_csv(snakemake.input.maf1, sep="\t", comment="#", low_memory=False, dtype=str)
maf2 = pd.read_csv(snakemake.input.maf2, sep="\t", comment="#", low_memory=False, dtype=str)

print(f"{caller1_name}: {len(maf1)} variants", file=log)
print(f"{caller2_name}: {len(maf2)} variants", file=log)

missing1 = [c for c in key_cols if c not in maf1.columns]
missing2 = [c for c in key_cols if c not in maf2.columns]
if missing1:
    raise ValueError(f"Key columns missing from {caller1_name} MAF: {missing1}")
if missing2:
    raise ValueError(f"Key columns missing from {caller2_name} MAF: {missing2}")

maf1["_key"] = maf1[key_cols].fillna("").agg("_".join, axis=1)
maf2["_key"] = maf2[key_cols].fillna("").agg("_".join, axis=1)

keys1 = set(maf1["_key"])
keys2 = set(maf2["_key"])

shared = keys1 & keys2
only1 = keys1 - keys2
only2 = keys2 - keys1

print(f"shared: {len(shared)}", file=log)
print(f"{caller1_name}-only: {len(only1)}", file=log)
print(f"{caller2_name}-only: {len(only2)}", file=log)

caller1_only_maf = maf1[maf1["_key"].isin(only1)].drop(columns=["_key"])
caller2_only_maf = maf2[maf2["_key"].isin(only2)].drop(columns=["_key"])
shared_maf = maf1[maf1["_key"].isin(shared)].drop(columns=["_key"])

caller1_only_maf.to_csv(snakemake.output.caller1_only, sep="\t", index=False, na_rep="")
caller2_only_maf.to_csv(snakemake.output.caller2_only, sep="\t", index=False, na_rep="")

# Driver gene analysis
driver_counts = {"intersect": 0, caller1_name: 0, caller2_name: 0}
driver_genes_found = {"intersect": "", caller1_name: "", caller2_name: ""}

if driver_genes_path:
    print(f"Loading driver gene list from {driver_genes_path}", file=log)
    driver_df = pd.read_csv(driver_genes_path, sep="\t", dtype=str)
    drivers = set(driver_df.loc[driver_df[driver_col] == driver_col_value, driver_gene_col])
    print(f"Loaded {len(drivers)} bona fide driver genes", file=log)

    for label, maf in [("intersect", shared_maf), (caller1_name, caller1_only_maf), (caller2_name, caller2_only_maf)]:
        if maf_gene_col not in maf.columns:
            print(f"WARNING: {maf_gene_col} not found in {label} MAF; skipping driver gene annotation", file=log)
            continue
        if maf_vc_col in maf.columns:
            coding = maf[maf[maf_vc_col].isin(coding_variant_classes)]
        else:
            print(f"WARNING: {maf_vc_col} not found in {label} MAF; counting all variant classes for driver genes", file=log)
            coding = maf
        hit_genes = set(coding[maf_gene_col].dropna()) & drivers
        hit_count = int(coding[maf_gene_col].isin(drivers).sum())
        genes_str = ",".join(sorted(hit_genes))
        suffix = "" if label == "intersect" else "-only"
        print(f"{label}{suffix} coding driver gene variants: {hit_count} across {len(hit_genes)} genes: {genes_str}", file=log)
        driver_counts[label] = hit_count
        driver_genes_found[label] = genes_str

# Long-format stats: one row per category
meta = {
    "tumour_id": snakemake.wildcards.tumour_id,
    "normal_id": snakemake.wildcards.normal_id,
    "seq_type": snakemake.wildcards.seq_type,
    "genome_build": snakemake.wildcards.genome_build,
    "pair_status": snakemake.wildcards.pair_status,
    "caller1": caller1_name,
    "caller2": caller2_name,
}
rows = [
    {**meta, "category": "intersect",
     "count": len(shared),
     "driver_variants": driver_counts["intersect"],
     "driver_genes": driver_genes_found["intersect"]},
    {**meta, "category": f"{caller1_name}_only",
     "count": len(only1),
     "driver_variants": driver_counts[caller1_name],
     "driver_genes": driver_genes_found[caller1_name]},
    {**meta, "category": f"{caller2_name}_only",
     "count": len(only2),
     "driver_variants": driver_counts[caller2_name],
     "driver_genes": driver_genes_found[caller2_name]},
]
stats = pd.DataFrame(rows)
stats.to_csv(snakemake.output.stats, sep="\t", index=False, na_rep="")

log.close()
