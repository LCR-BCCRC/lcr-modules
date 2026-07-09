import sys
import pandas as pd

log = open(snakemake.log[0], "w")
sys.stderr = log

caller1_name = snakemake.params.caller1_name
caller2_name = snakemake.params.caller2_name
key_cols = snakemake.params.key_cols

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

caller1_only_maf.to_csv(snakemake.output.caller1_only, sep="\t", index=False, na_rep="")
caller2_only_maf.to_csv(snakemake.output.caller2_only, sep="\t", index=False, na_rep="")

stats = pd.DataFrame([{
    "tumour_id": snakemake.wildcards.tumour_id,
    "normal_id": snakemake.wildcards.normal_id,
    "seq_type": snakemake.wildcards.seq_type,
    "genome_build": snakemake.wildcards.genome_build,
    "pair_status": snakemake.wildcards.pair_status,
    "caller1": caller1_name,
    "caller2": caller2_name,
    "caller1_total": len(maf1),
    "caller2_total": len(maf2),
    "shared": len(shared),
    "caller1_only": len(only1),
    "caller2_only": len(only2),
}])
stats.to_csv(snakemake.output.stats, sep="\t", index=False)

log.close()
