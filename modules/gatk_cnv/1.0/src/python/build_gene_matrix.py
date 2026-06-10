#!/usr/bin/env python3
"""Build per-gene relative and absolute CN matrices from per-tumour
GATK CallCopyRatioSegments output and PureCN absolute gene CN tables.

For each tumour:
  - Relative CN: from .called.seg, mapped to genes via a GTF (gene body overlap with
    segment, weighted by overlap length, then mean log2 ratio across the gene).
    Output value: 2^(weighted mean log2 ratio) — i.e. linear copy ratio relative to PoN baseline.
  - Absolute CN: prefer PureCN's <sample>_genes.csv (column "C") when populated; otherwise
    compute per-segment absolute CN from purity/ploidy + segment log2 ratio using the
    PureCN formula t = (2^L * (ρ*P + 2*(1-ρ)) - 2*(1-ρ)) / ρ, then take the
    overlap-length-weighted mean per gene and round to integer.

The fallback path is what gets used in this project: PureCN's callAlterations needs
a Gene column in the interval file we don't currently populate, so it returns no
per-gene calls. Computing absolute CN from purity/ploidy + segments matches what
callAlterations does internally.

Output: two TSVs, rows = genes, columns = sample_ids.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import Iterable

import numpy as np
import pandas as pd
import pyranges as pr


def load_gtf_genes(gtf_path: str) -> pd.DataFrame:
    """Parse gene records from a GTF, returning chrom/start/end/gene_name/gene_id."""
    cols = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    rows = []
    open_fn = open
    if gtf_path.endswith(".gz"):
        import gzip
        open_fn = gzip.open
    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = {}
            for kv in parts[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                k, _, v = kv.partition(" ")
                attrs[k] = v.strip().strip('"')
            rows.append({
                "chrom": parts[0].lstrip("chr"),
                "start": int(parts[3]) - 1,  # GTF is 1-based inclusive; convert to 0-based
                "end":   int(parts[4]),
                "gene_id":   attrs.get("gene_id", ""),
                "gene_name": attrs.get("gene_name", attrs.get("gene_id", "")),
            })
    df = pd.DataFrame(rows)
    df = df.dropna(subset=["gene_name"]).drop_duplicates(subset=["gene_name"])
    return df


def load_called_seg(seg_path: str) -> pd.DataFrame:
    """Parse a GATK CallCopyRatioSegments .called.seg file.
    Columns we care about: CONTIG, START, END, NUM_POINTS_COPY_RATIO, MEAN_LOG2_COPY_RATIO, CALL.
    Skip @ header lines."""
    rows = []
    with open(seg_path) as fh:
        header = None
        for line in fh:
            if line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            try:
                rows.append({
                    "chrom": row["CONTIG"].lstrip("chr"),
                    "start": int(row["START"]) - 1,
                    "end":   int(row["END"]),
                    "log2": float(row["MEAN_LOG2_COPY_RATIO"]),
                })
            except (KeyError, ValueError):
                continue
    return pd.DataFrame(rows)


def _join_genes_segs(genes: pd.DataFrame, segs: pd.DataFrame) -> pd.DataFrame:
    """pyranges-vectorized gene × segment overlap. Returns one row per overlap
    with gene_name + segment columns + an `overlap_len` column (in bp)."""
    g_pr = pr.PyRanges(genes.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
    s_pr = pr.PyRanges(segs.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
    j = g_pr.join(s_pr).df
    if j.empty:
        return j
    j["overlap_len"] = (j[["End", "End_b"]].min(axis=1) - j[["Start", "Start_b"]].max(axis=1)).clip(lower=0)
    return j


def gene_log2_from_segments(genes: pd.DataFrame, segs: pd.DataFrame) -> pd.Series:
    """For each gene, compute overlap-length-weighted mean log2 across overlapping segments.
    Returns a Series indexed by gene_name. Vectorized via pyranges."""
    if segs.empty:
        return pd.Series(np.nan, index=genes["gene_name"].values, name="log2")
    j = _join_genes_segs(genes, segs)
    if j.empty:
        return pd.Series(np.nan, index=genes["gene_name"].values, name="log2")
    j["w"] = j["overlap_len"] * j["log2"]
    agg = j.groupby("gene_name").agg(w=("w", "sum"), L=("overlap_len", "sum"))
    out = (agg["w"] / agg["L"]).where(agg["L"] > 0, other=np.nan)
    return out.reindex(genes["gene_name"].values).rename("log2")


def sample_id_from_path(path: str, suffix: str) -> str:
    base = os.path.basename(path)
    if base.endswith(suffix):
        base = base[: -len(suffix)]
    # Path like .../<seq_type>--<genome_build>/<tumour_id>--<normal_id>--<pair_status>.called.seg
    # The leading "<tumour_id>" is what we want as the column name.
    return base.split("--")[0]


def load_purecn_genes(genes_csv: str) -> pd.Series:
    """Load PureCN's per-gene absolute CN table. Returns Series gene_name -> integer CN.
    Empty file or missing 'C' column → empty Series (caller falls back to compute_absolute_cn)."""
    if not os.path.exists(genes_csv) or os.path.getsize(genes_csv) == 0:
        return pd.Series(dtype=float, name="C")
    df = pd.read_csv(genes_csv)
    # PureCN's callAlterations() output indexes by gene name in the row name.
    # When written via write.csv with row.names=TRUE, the first column is unnamed — pandas reads it as "Unnamed: 0".
    if "gene.symbol" in df.columns:
        gene_col = "gene.symbol"
    elif "Unnamed: 0" in df.columns:
        gene_col = "Unnamed: 0"
    elif "gene_name" in df.columns:
        gene_col = "gene_name"
    else:
        gene_col = df.columns[0]
    if "C" not in df.columns:
        return pd.Series(dtype=float, name="C")
    return pd.Series(df["C"].values, index=df[gene_col].astype(str), name="C")


def load_purity_csv(purity_csv: str) -> tuple[float, float]:
    """Read purity and ploidy from a PureCN purity CSV (one row per sample).
    Returns (purity, ploidy); (NaN, NaN) if file is missing/empty/malformed."""
    if not os.path.exists(purity_csv) or os.path.getsize(purity_csv) == 0:
        return float("nan"), float("nan")
    try:
        df = pd.read_csv(purity_csv)
    except Exception:
        return float("nan"), float("nan")
    if df.empty or "purity" not in df.columns or "ploidy" not in df.columns:
        return float("nan"), float("nan")
    return float(df["purity"].iloc[0]), float(df["ploidy"].iloc[0])


def absolute_cn_from_log2(log2_ratio: float, purity: float, ploidy: float) -> float:
    """Continuous absolute tumor copy number from segment log2 ratio + purity/ploidy.
    PureCN formula: t = (2^L * (ρ*P + 2*(1-ρ)) - 2*(1-ρ)) / ρ.
    Returns NaN if inputs are invalid; clamps result to >= 0."""
    if pd.isna(log2_ratio) or pd.isna(purity) or pd.isna(ploidy) or purity <= 0:
        return float("nan")
    cr = 2.0 ** log2_ratio
    t = (cr * (purity * ploidy + 2 * (1 - purity)) - 2 * (1 - purity)) / purity
    return max(0.0, float(t))


def gene_absolute_cn_from_segments(genes: pd.DataFrame, segs: pd.DataFrame,
                                   purity: float, ploidy: float) -> pd.Series:
    """Per-gene integer absolute CN: weighted mean of per-segment continuous CN over
    overlapping segments, then rounded. Returns Series gene_name -> int CN (NaN if
    no overlapping segments or purity/ploidy unavailable). Vectorized via pyranges.

    Self-consistency centering: PureCN's reported purity/ploidy is fit against
    re-centered log2 values (we typically see offsets of 0.05-0.1 log2). If we
    apply the formula directly to GATK's raw log2, the implied sample-average CN
    won't equal the reported ploidy. We compute the length-weighted mean log2
    across the sample and subtract it from each segment, so by construction
    the weighted-mean integer CN matches `ploidy`."""
    if segs.empty or pd.isna(purity) or pd.isna(ploidy) or purity <= 0:
        return pd.Series(np.nan, index=genes["gene_name"].values, name="C")

    seg_lens_all = (segs["end"] - segs["start"]).clip(lower=0)
    total_len = float(seg_lens_all.sum())
    offset = float((segs["log2"] * seg_lens_all).sum() / total_len) if total_len > 0 else 0.0

    # Compute per-segment continuous CN once (vectorized).
    cr = 2.0 ** (segs["log2"] - offset)
    seg_cn = (cr * (purity * ploidy + 2 * (1 - purity)) - 2 * (1 - purity)) / purity
    seg_cn = seg_cn.clip(lower=0)

    segs_with_cn = segs.copy()
    segs_with_cn["seg_cn"] = seg_cn

    j = _join_genes_segs(genes, segs_with_cn)
    if j.empty:
        return pd.Series(np.nan, index=genes["gene_name"].values, name="C")
    j["w"] = j["overlap_len"] * j["seg_cn"]
    agg = j.groupby("gene_name").agg(w=("w", "sum"), L=("overlap_len", "sum"))
    out = (agg["w"] / agg["L"]).where(agg["L"] > 0, other=np.nan).round()
    return out.reindex(genes["gene_name"].values).rename("C")


def build_matrix(per_sample: dict, all_genes: list) -> pd.DataFrame:
    return pd.DataFrame(per_sample, index=all_genes)


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True, help="Ensembl/Gencode GTF for the matching build")
    p.add_argument("--called-segs", required=True, nargs="+", help="GATK .called.seg files (one per tumour)")
    p.add_argument("--purecn-genes", required=True, nargs="+", help="PureCN <sample>_genes.csv files (one per tumour)")
    p.add_argument("--purity-csvs", required=True, nargs="+", help="PureCN <sample>.csv purity/ploidy files (one per tumour)")
    p.add_argument("--relative-out", required=True, help="Output TSV: gene × sample relative CN (linear ratio)")
    p.add_argument("--absolute-out", required=True, help="Output TSV: gene × sample absolute CN (integer)")
    args = p.parse_args(list(argv) if argv is not None else None)

    print(f"Loading GTF: {args.gtf}", file=sys.stderr)
    genes = load_gtf_genes(args.gtf)
    if genes.empty:
        print("ERROR: no gene records found in GTF", file=sys.stderr)
        return 1
    print(f"  → {len(genes)} unique gene records", file=sys.stderr)

    # Build per-sample lookups for purity/ploidy and PureCN per-gene calls.
    purity_by_sid = {}
    for pc in args.purity_csvs:
        sid = sample_id_from_path(pc, ".csv")
        purity_by_sid[sid] = load_purity_csv(pc)
    purecn_genes_by_sid = {sample_id_from_path(g, "_genes.csv"): g for g in args.purecn_genes}

    relative_per_sample = {}
    absolute_per_sample = {}
    for seg_path in args.called_segs:
        sid = sample_id_from_path(seg_path, ".called.seg")
        segs = load_called_seg(seg_path)

        # Relative CN: weighted mean log2 → linear copy ratio.
        rel_log2 = gene_log2_from_segments(genes, segs)
        relative_per_sample[sid] = (2.0 ** rel_log2).reindex(genes["gene_name"].values)

        # Absolute CN: prefer PureCN's per-gene calls when populated; otherwise compute
        # from segment log2 + purity/ploidy.
        purecn_genes_path = purecn_genes_by_sid.get(sid)
        purecn_series = load_purecn_genes(purecn_genes_path) if purecn_genes_path else pd.Series(dtype=float, name="C")
        if not purecn_series.empty:
            print(f"  + absolute CN for {sid}: PureCN per-gene calls ({len(purecn_series)} genes)", file=sys.stderr)
            absolute_per_sample[sid] = purecn_series.reindex(genes["gene_name"].values)
        else:
            purity, ploidy = purity_by_sid.get(sid, (float("nan"), float("nan")))
            print(f"  + absolute CN for {sid}: computed from segments (purity={purity}, ploidy={ploidy})", file=sys.stderr)
            abs_series = gene_absolute_cn_from_segments(genes, segs, purity, ploidy)
            absolute_per_sample[sid] = abs_series.reindex(genes["gene_name"].values)

    rel_df = pd.DataFrame(relative_per_sample, index=genes["gene_name"].values)
    abs_df = pd.DataFrame(absolute_per_sample, index=genes["gene_name"].values)

    rel_df.index.name = "gene_name"
    abs_df.index.name = "gene_name"

    rel_df.to_csv(args.relative_out, sep="\t", na_rep="NA")
    abs_df.to_csv(args.absolute_out, sep="\t", na_rep="NA")

    print(f"Wrote {args.relative_out}: {rel_df.shape[0]} genes × {rel_df.shape[1]} samples", file=sys.stderr)
    print(f"Wrote {args.absolute_out}: {abs_df.shape[0]} genes × {abs_df.shape[1]} samples", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
