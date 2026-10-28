import sys
import pandas as pd

log = open(snakemake.log[0], "w")
sys.stderr = log

driver_genes_path = snakemake.params.driver_genes
driver_gene_col = snakemake.params.driver_gene_col
driver_col = snakemake.params.driver_col
driver_col_value = snakemake.params.driver_col_value
maf_gene_col = snakemake.params.maf_gene_col
maf_vc_col = snakemake.params.maf_vc_col
coding_variant_classes = set(snakemake.params.coding_variant_classes)

maf = pd.read_csv(snakemake.input.maf, sep="\t", comment="#", low_memory=False, dtype=str)
print(f"Input MAF: {len(maf)} variants", file=log)

# Load driver gene list
driver_df = pd.read_csv(driver_genes_path, sep="\t", dtype=str)
drivers = set(driver_df.loc[driver_df[driver_col] == driver_col_value, driver_gene_col])
print(f"Loaded {len(drivers)} bona fide driver genes", file=log)

# Filter to coding variants in driver genes
if maf_vc_col not in maf.columns:
    raise ValueError(f"Variant classification column '{maf_vc_col}' not found in MAF")
if maf_gene_col not in maf.columns:
    raise ValueError(f"Gene column '{maf_gene_col}' not found in MAF")

coding = maf[maf[maf_vc_col].isin(coding_variant_classes)]
filtered = coding[coding[maf_gene_col].isin(drivers)]

print(f"After coding filter: {len(coding)} variants", file=log)
print(f"After driver gene filter: {len(filtered)} variants", file=log)
genes_hit = sorted(filtered[maf_gene_col].dropna().unique())
print(f"Driver genes represented: {','.join(genes_hit)}", file=log)

info_columns = snakemake.params.info_columns
missing = [c for c in info_columns if c not in maf.columns]
if missing:
    print(f"WARNING: the following info_columns are not present in the MAF and will be ignored by igv-reports: {missing}", file=log)

# Tabulator.js uses regex matching on cell values; unescaped regex metacharacters
# (e.g. '?' in HGVS notation like 'p.M1?') crash table initialization silently.
# Replace them in display columns only — the rest of the MAF is unaffected.
display_cols = [c for c in info_columns if c in filtered.columns and filtered[c].dtype == object]
filtered = filtered.copy()
for col in display_cols:
    filtered[col] = filtered[col].str.replace("?", "_loss", regex=False)

filtered.to_csv(snakemake.output.maf, sep="\t", index=False, na_rep="")

log.close()
