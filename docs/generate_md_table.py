import pandas as pd
import plotly.graph_objects as go
import re

input_file = "module_seq_types.tsv"
svg_output = "modules_alluvial.svg"

base_url = "https://github.com/LCR-BCCRC/lcr-modules/tree/master/modules/"
readme_url = "https://github.com/LCR-BCCRC/lcr-modules/tree/demo_repairs?tab=readme-ov-file"

def github_anchor(text):
    text = text.lower()
    text = re.sub(r"[^\w\s-]", "", text)
    text = re.sub(r"\s+", "-", text)
    return text


modules = pd.read_csv(input_file, sep="\t")

modules["purpose_clean"] = (
    modules["purpose"]
    .str.replace(r"\s*\(N=\d+\)$", "", regex=True)
)

short_reads = ["genome", "capture", "mrna"]

modules["data_type"] = "Long reads"
modules.loc[
    modules["seq_type"].isin(short_reads),
    "data_type"
] = "Illumina short reads"

counts = (
    modules
    .groupby(["data_type", "seq_type", "purpose_clean"])
    .size()
    .reset_index(name="value")
)

# counts per purpose for labels
purpose_counts = (
    modules
    .groupby("purpose_clean")
    .size()
    .to_dict()
)

data_nodes = sorted(counts["data_type"].unique())
seq_nodes = sorted(counts["seq_type"].unique())
purpose_nodes = sorted(counts["purpose_clean"].unique())

nodes = data_nodes + seq_nodes + purpose_nodes
node_index = {n: i for i, n in enumerate(nodes)}

labels = []

for n in nodes:

    if n in purpose_nodes:

        count = purpose_counts[n]
        anchor = github_anchor(n)

        labels.append(
            f"<a href='{readme_url}#{anchor}'>{n} (N={count})</a>"
        )

    else:
        labels.append(n)

sources = []
targets = []
values = []

# data_type -> seq_type
d1 = (
    modules
    .groupby(["data_type", "seq_type"])
    .size()
    .reset_index(name="value")
)

for _, r in d1.iterrows():

    sources.append(node_index[r["data_type"]])
    targets.append(node_index[r["seq_type"]])
    values.append(r["value"])

# seq_type -> purpose
d2 = (
    modules
    .groupby(["seq_type", "purpose_clean"])
    .size()
    .reset_index(name="value")
)

for _, r in d2.iterrows():

    sources.append(node_index[r["seq_type"]])
    targets.append(node_index[r["purpose_clean"]])
    values.append(r["value"])

color_map = {
    "Illumina short reads": "#ff7f0e",
    "Long reads": "#9467bd",
    "genome": "#1f77b4",
    "capture": "#2ca02c",
    "mrna": "#d62728",
    "promethION": "#17becf"
}

node_colors = []

for n in nodes:

    if n in color_map:
        node_colors.append(color_map[n])
    else:
        node_colors.append("#bbbbbb")

fig = go.Figure(
    go.Sankey(
        arrangement="snap",

        node=dict(
            pad=25,
            thickness=25,
            line=dict(color="black", width=0.3),
            label=labels,
            color=node_colors
        ),

        link=dict(
            source=sources,
            target=targets,
            value=values,
            color="rgba(150,150,150,0.35)"
        )
    )
)

fig.update_layout(
    font=dict(size=14),
    margin=dict(l=20, r=20, t=20, b=20),
    width=1200,
    height=900
)

fig.write_image(svg_output)

## Generate markdown for readme

# Convert module names to markdown links
modules["module"] = (
    modules["module"]
    .apply(lambda m: f"[{m}]({base_url}{m})")
)

# Columns to include in tables
table_cols = [c for c in modules.columns if c not in ["purpose", "purpose_clean"]]

modules = (
    modules
    .groupby(["purpose_clean", "module"], as_index=False)
    .agg({
        "seq_type": lambda x: "; ".join(sorted(x)),
        "data_type": "first"
    })
)

def github_anchor(text):
    text = text.lower()
    text = re.sub(r"[^\w\s-]", "", text)
    text = re.sub(r"\s+", "-", text)
    return text

with open(output_file, "w") as f:

    # Write TOC
    for p in sorted(modules["purpose_clean"].unique()):
        f.write(f"- [{p}](#{github_anchor(p)})\n")

    f.write("\n---\n\n")

    for purpose, group in modules.groupby("purpose_clean"):
        f.write(f"### {purpose}\n\n")

        # Table header
        header = "| " + " | ".join(table_cols) + " |\n"
        sep = "|" + "|".join(["---"] * len(table_cols)) + "|\n"

        f.write(header)
        f.write(sep)

        # Table rows
        for _, row in group.iterrows():
            row_vals = [str(row[c]) for c in table_cols]
            f.write("| " + " | ".join(row_vals) + " |\n")

        f.write("\n")
