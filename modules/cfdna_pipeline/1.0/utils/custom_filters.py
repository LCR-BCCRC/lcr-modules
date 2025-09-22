"""
Variants are filtered through multiple steps:

1) Allele frequency is recalculated based on t_alt_count/t_depth.
2) gnomAD_AF is filtered by a threshold across all population columns.
3) Variants are marked as potential CHIP if they have VAF in normal >2% or
   normal VAF is >2x (or >3x for non-CHIP genes) of tumor VAF.
4) Variants are filtered out if they have less than min_alt_tum in the tumour
   or less than min_germline_depth in the germline.
5) Additional filter removes variants with low normal depth (<=100) and high AF (>0.3).
6) Minimum tumor depth filter is applied.
7) Variants are marked if they are in a blacklisted position or a hotspot.
   Blacklisted variants are removed.
8) Variants are filtered out if they have VAF < min_tumour_vaf unless they are
   part of a phase set (>1 variant) or in a hotspot.
9) Minimum UMI support filter is applied (tiered thresholds on UMI_3_count).
10) A variant_source column is added indicating: phase_group, high_vaf, hotspot,
    or combinations thereof.
"""
from dataclasses import dataclass
from collections import Counter
from typing import Set
import argparse
import pandas as pd

GNOMAD_COLS = [
    "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
    "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
]
BAD_TOKENS = {"", "na", "n/a", "none", "nan"}

def get_args():
    """Define and parse CLI arguments for the variant filtering pipeline."""
    p = argparse.ArgumentParser(
        description="Filter variants through gnomAD, CHIP, blacklist/hotspots, VAF/phase, and read support rules."
    )
    # I/O
    p.add_argument('--input_maf', required=True, type=str, help='Path to input MAF file (tab-delimited).')
    p.add_argument('--output_maf', required=True, type=str, help='Path to write the filtered MAF.')
    p.add_argument('--blacklist', required=True, type=str, help='TSV with blacklisted variant_key (Chromosome:Start_Position).')
    p.add_argument('--hotspots', required=True, type=str, help='TSV with hotspot variant_key (Chromosome:Start_Position).')
    # thresholds
    p.add_argument('--gnomad_threshold', required=True, type=float, help='Max allowed gnomAD AF across all population columns.')
    p.add_argument('--min_alt_depth_tum', required=True, type=int, help='Minimum tumor alt read count (t_alt_count).')
    p.add_argument('--min_germline_depth', required=True, type=int, help='Minimum normal depth (n_depth).')
    p.add_argument('--min_tumour_vaf', type=float, default=0.01, help='Minimum VAF to keep unless in phase group or hotspot.')
    p.add_argument('--min_t_depth', type=int, default=200, help='Minimum tumor depth (t_depth).')
    # UMI thresholds
    p.add_argument('--min_UMI_3_count', type=int, default=1, help='Minimum UMI_3_count for variants at/above low_alt_thresh.')
    p.add_argument('--low_alt_thresh', type=int, default=20, help='Alt read threshold below which stricter UMI_3_count applies.')
    p.add_argument('--low_alt_min_UMI_3_count', type=int, default=2, help='Minimum UMI_3_count for low-alt variants (< low_alt_thresh).')
    # low-normal safeguards
    p.add_argument('--low_normal_depth', type=int, default=100, help='If n_depth <= this, apply stricter AF filter.')
    p.add_argument('--low_normal_AF', type=float, default=0.30, help='AF threshold used when n_depth is low.')
    return p.parse_args()

@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for the variant filtering pipeline."""
    # I/O
    input_maf: str
    output_maf: str
    blacklist: str
    hotspots: str
    # thresholds
    gnomad_threshold: float
    min_alt_depth_tum: int
    min_germline_depth: int
    min_tumour_vaf: float = 0.01
    min_t_depth: int = 200
    # UMI thresholds
    min_UMI_3_count: int = 1
    low_alt_thresh: int = 20
    low_alt_min_UMI_3_count: int = 2
    # low-normal safeguards
    low_normal_depth: int = 100
    low_normal_AF: float = 0.30

class VariantFilterPipeline:
    """Variant filtering pipeline object."""
    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg

    def filter_variants(self) -> pd.DataFrame:
        """Run the full pipeline in order and return the filtered variants."""
        df = self.read_maf(self.cfg.input_maf)
        # if no variants skip rest of pipeline
        if df.empty:
            return df
        df = self.recalculate_af(df)
        df = self.filter_gnomad(df, self.cfg.gnomad_threshold)
        df = self.mark_potential_chip(df)
        df = self.mark_blacklist_hotspot(df, self.cfg.blacklist, self.cfg.hotspots)
        df = self.remove_blacklisted(df)
        df = self.filter_vaf_and_phase(df, self.cfg.min_tumour_vaf)
        df = self.filter_by_read_support(
            df,
            min_t_depth=self.cfg.min_t_depth,
            min_germline_depth=self.cfg.min_germline_depth,
            min_alt_tum=self.cfg.min_alt_depth_tum,
            min_UMI_3_count=self.cfg.min_UMI_3_count,
            low_alt_thresh=self.cfg.low_alt_thresh,
            low_alt_UMI_3_count=self.cfg.low_alt_min_UMI_3_count,
            low_normal_depth=self.cfg.low_normal_depth,
            low_normal_AF=self.cfg.low_normal_AF,
        )
        return df

    # Steps (each takes df -> returns df)
    def read_maf(self, maf_file: str) -> pd.DataFrame:
        """Read a MAF file and add a variant_key column (Chromosome:Start_Position)."""
        df = pd.read_csv(maf_file, sep="\t")
        df["variant_key"] = df["Chromosome"].astype(str) + ":" + df["Start_Position"].astype(str)
        return df

    def recalculate_af(self, df: pd.DataFrame) -> pd.DataFrame:
        """Recalculate AF as t_alt_count / t_depth."""
        out = df.copy()
        out["AF"] = out["t_alt_count"] / out["t_depth"]
        return out

    def filter_gnomad(self, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        """Keep rows where all available gnomAD population AF columns are < threshold or NaN."""
        cols = [c for c in GNOMAD_COLS if c in df.columns]
        if not cols:
            return df.copy()
        keep = ((df[cols] < threshold) | df[cols].isna()).all(axis=1)
        return df.loc[keep].copy()

    def mark_blacklist_hotspot(self, df: pd.DataFrame, blacklist: str, hotspots: str) -> pd.DataFrame:
        """Add blacklist and hotspot boolean flags based on variant_key membership."""
        out = df.copy()
        bl = pd.read_csv(blacklist, sep="\t"); bl.columns = ["variant_key"]
        hs = pd.read_csv(hotspots, sep="\t"); hs.columns = ["variant_key"]
        out["blacklist"] = out["variant_key"].isin(bl["variant_key"])
        out["hotspot"] = out["variant_key"].isin(hs["variant_key"])
        return out

    def remove_blacklisted(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drop variants flagged as blacklisted (if flag is present)."""
        if "blacklist" not in df.columns:
            return df.copy()
        return df.loc[df["blacklist"] == False].copy()

    def _fetch_phase_sets(self, lps_series: pd.Series) -> Set[str]:
        """Return LPS IDs that occur in >1 variant (handles comma-separated multi-membership)."""
        s = lps_series.fillna("").astype(str).str.split(",").explode().str.strip()
        s = s[~s.str.lower().isin(BAD_TOKENS)]
        counts = Counter(s)
        return {ps for ps, c in counts.items() if c > 1}

    def filter_vaf_and_phase(self, df: pd.DataFrame, min_vaf: float) -> pd.DataFrame:
        """Keep variants if in a multi-member phase set, have AF >= min_vaf, or are hotspots. Label variant_source."""
        out = df.copy()
        out["variant_source"] = [[] for _ in range(out.shape[0])]
        phase_sets = self._fetch_phase_sets(out["LPS"]) if "LPS" in out.columns else set()

        for row in out.itertuples():
            lsps = [t.strip() for t in str(getattr(row, "LPS", "")).split(",")]
            lsps = [t for t in lsps if t and t.lower() not in BAD_TOKENS]
            if phase_sets and any(t in phase_sets for t in lsps):
                row.variant_source.append("phase_group")
            if getattr(row, "AF", 0) >= min_vaf:
                row.variant_source.append("high_vaf")
            if bool(getattr(row, "hotspot", False)):
                row.variant_source.append("hotspot")

        out["variant_source"] = out["variant_source"].apply(lambda xs: ",".join(xs) if xs else None)
        out = out.loc[out["variant_source"].notna()].reset_index(drop=True).copy()
        return out

    def mark_potential_chip(self, df: pd.DataFrame) -> pd.DataFrame:
        """Mark likely CHIP variants based on VAF in normal vs tumor and stricter rules for CHIP genes."""
        chip_genes = {"DNMT3A","TET2","ASXL1","PPM1D","TP53","JAK2","SF3B1","SRSF2"}
        out = df.copy()
        out["VAF_normal"] = (out["n_alt_count"] / out["n_depth"].replace(0, pd.NA)).fillna(0.0)
        AF = out.get("AF", (out["t_alt_count"] / out["t_depth"].replace(0, pd.NA))).fillna(0.0)
        normal_gt_tumor = (out["VAF_normal"] / AF) >= 3
        normal_high = out["VAF_normal"] >= 0.02
        chip_like = normal_gt_tumor | normal_high
        chip_strict = chip_like | (out["Hugo_Symbol"].isin(chip_genes) & (out["n_alt_count"] >= 1))
        out["CHIP"] = chip_strict
        return out

    def filter_by_read_support(
        self,
        df: pd.DataFrame,
        min_t_depth: int,
        min_germline_depth: int,
        min_alt_tum: int,
        min_UMI_3_count: int = 1,
        low_alt_thresh: int = 20,
        low_alt_UMI_3_count: int = 2,
        low_normal_depth: int = 100,
        low_normal_AF: float = 0.30,
    ) -> pd.DataFrame:
        """Apply tumor/normal depth, AF under low-normal-depth, and tiered UMI_3_count filters; always keep hotspots."""
        cols = {"t_depth","n_depth","t_alt_count","UMI_3_count","AF","hotspot"}
        missing = cols - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

        mask = (
            (df["t_alt_count"] >= min_alt_tum)
            & (df["n_depth"] >= min_germline_depth)
            & (df["t_depth"] >= min_t_depth)
        )
        mask &= ~((df["n_depth"] <= low_normal_depth) & (df["AF"] > low_normal_AF))
        umi_mask = (
            ((df["t_alt_count"] >= low_alt_thresh) & (df["UMI_3_count"] >= min_UMI_3_count))
            | ((df["t_alt_count"] < low_alt_thresh) & (df["UMI_3_count"] >= low_alt_UMI_3_count))
        )
        mask &= umi_mask
        mask |= df["hotspot"].eq(True)
        return df.loc[mask].copy()

def main():
    """Parse args, run the pipeline, and write the filtered MAF."""
    args = get_args()
    cfg = PipelineConfig(**vars(args))
    var_object = VariantFilterPipeline(cfg)

    df = var_object.filter_variants()
    print(f"Finished custom filtering, {df.shape[0]} variants remain")
    df.to_csv(cfg.output_maf, sep="\t", index=False)

if __name__ == '__main__':
    main()
