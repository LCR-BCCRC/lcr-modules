"""
Variants are filtered through multiple steps:

1) Allele frequency (AF) is recalculated based on t_alt_count / t_depth.
2) Variants are marked as potential CHIP (annotation only).
3) Blacklist and hotspot positions annotated.
4) Phasing annotated from LPS: is_phased (>=2 members in a phase group), phase_set_size.
   hotspot or phased variants are considered "exempt" from subsequent threshold-based filtering.
5) gnomAD filtering: keep if all gnomAD AF population columns < threshold OR hotspot OR is_phased.
6) Blacklisted variants are removed (blacklist always overrides exemptions).
7) VAF filtering: only non-hotspot, non-phased variants are removed if AF < min_tumour_vaf.
   All variants meeting the AF threshold (including hotspot/phased) get the "high_vaf" tag.
8) Read support filtering (tumor alt count, tumor depth, normal depth, low-normal AF guard, tiered UMI thresholds).
   hotspot OR phased variants bypass these filters entirely.
9) variant_source accumulates all applicable reasons: hotspot, phase_group, high_vaf.
10) Final filtered MAF written.

Exemption rule summary:
A variant that is hotspot OR phased is retained regardless of gnomAD, VAF, read support, depth, or UMI thresholds (unless blacklisted).
"""
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import argparse
import pandas as pd
from collections import Counter
import os
import gzip

GNOMAD_COLS = [
    "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
    "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
]
BAD_TOKENS = {"", "na", "n/a", "none", "nan"}

def get_args():
    """Define and parse CLI arguments for the variant filtering pipeline."""
    p = argparse.ArgumentParser(
        description="Filter variants with hotspot/phasing exemptions and cumulative variant_source tagging."
    )
    # I/O
    p.add_argument('--input_maf', required=True, type=str, help='Path to input MAF file (tab-delimited).')
    p.add_argument('--output_maf', required=True, type=str, help='Path to write the filtered MAF.')
    p.add_argument('--blacklist', type=str, default="", help='TSV with blacklisted variant_key (Chromosome:Start_Position). Optional.')
    p.add_argument('--hotspots', type=str, default="", help='VCF file with hotspot variants. Optional.')
    # thresholds
    p.add_argument('--gnomad_threshold', required=False, type=float, help='Max allowed gnomAD AF across all population columns.')
    p.add_argument('--min_alt_depth_tum', required=False, type=int, help='Minimum tumor alt read count (t_alt_count).')
    p.add_argument('--min_germline_depth', required=False, type=int, help='Minimum normal depth (n_depth).')
    p.add_argument('--min_tumour_vaf', type=float, default=0.01, help='Minimum VAF to keep non-hotspot/non-phased variants.')
    p.add_argument('--min_t_depth', type=int, default=200, help='Minimum tumor depth (t_depth).')
    # UMI thresholds
    p.add_argument('--min_UMI_3_count', type=int, default=1, help='Minimum UMI_3_count for variants at/above low_alt_thresh.')
    p.add_argument('--low_alt_thresh', type=int, default=20, help='Alt read threshold below which stricter UMI_3_count applies.')
    p.add_argument('--low_alt_min_UMI_3_count', type=int, default=2, help='Minimum UMI_3_count for low-alt variants (< low_alt_thresh).')
    # low-normal safeguards
    p.add_argument('--low_normal_depth', type=int, default=100, help='If n_depth <= this, apply stricter AF filter (unless exempt).')
    p.add_argument('--low_normal_AF', type=float, default=0.30, help='AF threshold used when n_depth is low (unless exempt).')
    # backgrount mutation rate filter
    p.add_argument('--background_rates', type=str, default="", 
               help='Optional TSV with background mutation rates per position.')
    p.add_argument('--min_background_samples', type=int, default=20, 
                help='Minimum n_samples required in background index.')
    p.add_argument('--background_n_std', type=float, default=2.0, 
                help='Number of standard deviations above mean for background filter.')
    return p.parse_args()

@dataclass(frozen=True)
class PipelineConfig:
    """Configuration for the variant filtering pipeline."""
    input_maf: str
    output_maf: str
    blacklist: str = ""
    hotspots: str = ""
    gnomad_threshold: float = 0.01
    min_alt_depth_tum: int = 5
    min_germline_depth: int = 50
    min_tumour_vaf: float = 0.01
    min_t_depth: int = 200
    min_UMI_3_count: int = 1
    low_alt_thresh: int = 20
    low_alt_min_UMI_3_count: int = 2
    low_normal_depth: int = 100
    low_normal_AF: float = 0.30
    background_rates: str = ""
    min_background_samples: int = 20
    background_n_std: float = 2.0

def read_hotspot_vcf(vcf_file: str) -> set:
    """
    Read a VCF file and return a set of (chrom, pos, alt) tuples.
    Skips header lines starting with '#'.
    For multi-allelic sites (comma-separated ALT), creates separate entries for each ALT.
    Returns empty set if file path is empty or file doesn't exist.
    Handles both plain text and gzipped VCF files.
    """
    
    if not vcf_file or not os.path.isfile(vcf_file):
        return set()
    
    hotspot_set = set()
    
    # Detect if file is gzipped
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'  # text mode for gzip
    else:
        opener = open
        mode = 'r'
    
    with opener(vcf_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            pos = fields[1]
            alt_field = fields[4]  # ALT column
            
            # Handle multi-allelic sites
            for alt in alt_field.split(','):
                hotspot_set.add((chrom, pos, alt.strip()))
    
    return hotspot_set

class VariantFilterPipeline:
    """Variant filtering pipeline with hotspot/phasing exemptions and cumulative variant_source labeling."""

    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg

    def filter_variants(self) -> pd.DataFrame:
        """Run the full variant filtering pipeline."""
        df = self.read_maf(self.cfg.input_maf)
        if df.empty:
            return df
        df = self.recalculate_af(df)
        df = self.mark_potential_chip(df)
        df = self.mark_blacklist_hotspot(df, self.cfg.blacklist, self.cfg.hotspots)
        df = self.mark_phased_vars(df)  # adds is_phased, phase_set_size, seeds variant_source with hotspot/phase_group
        df = self.filter_gnomad(df, self.cfg.gnomad_threshold)
        df = self.remove_blacklisted(df)
        df = self.filter_vaf(df, self.cfg.min_tumour_vaf)
        df = self.filter_background_rates(
            df,
            self.cfg.background_rates,
            self.cfg.min_background_samples,
            self.cfg.background_n_std
        )
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
        # Finalize variant_source formatting (ensure comma-separated strings)
        df["variant_source"] = df["variant_source"].apply(
            lambda v: ",".join(v) if isinstance(v, list) else v
        )
        return df

    # --- Core steps ---
    def read_maf(self, maf_file: str) -> pd.DataFrame:
        """Read a MAF file and add a variant_key column (Chromosome:Start_Position)."""
        df = pd.read_csv(maf_file, sep="\t")
        df["variant_key"] = df["Chromosome"].astype(str) + ":" + df["Start_Position"].astype(str)
        # remove 0 count alleles that can be placeholders, like DNPs and TNPs
        df["t_alt_count"] = pd.to_numeric(df["t_alt_count"], errors="coerce")
        df = df.loc[df["t_alt_count"].gt(0)].reset_index(drop=True)

        return df.copy()

    def recalculate_af(self, df: pd.DataFrame) -> pd.DataFrame:
        """Recalculate AF as t_alt_count / t_depth."""
        out = df.copy()
        out["AF"] = out["t_alt_count"] / out["t_depth"]
        return out

    def mark_potential_chip(self, df: pd.DataFrame) -> pd.DataFrame:
        """Annotate CHIP-like variants based on normal vs tumor VAF and CHIP gene presence."""
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

    def mark_blacklist_hotspot(self, df: pd.DataFrame, blacklist: str, hotspots: str) -> pd.DataFrame:
        """
        Add blacklist and hotspot boolean flags.
        Blacklist: matched on variant_key (Chromosome:Start_Position). Optional.
        Hotspot: matched on (Chromosome, Start_Position, Tumor_Seq_Allele2) from VCF. Optional.
        If file paths are empty or files don't exist, all variants are marked False for that category.
        """
        out = df.copy()
        
        # Read blacklist (TSV with variant_key) - optional
        if blacklist and os.path.isfile(blacklist):
            bl = pd.read_csv(blacklist, sep="\t")
            bl.columns = ["variant_key"]
            out["blacklist"] = out["variant_key"].isin(bl["variant_key"])
        else:
            out["blacklist"] = False
        
        # Read hotspots from VCF - optional
        hotspot_set = read_hotspot_vcf(hotspots)
        
        if hotspot_set:
            # Create matching key for MAF: (Chromosome, Start_Position, Alt_allele)
            out["hotspot"] = out.apply(
                lambda row: (str(row["Chromosome"]), str(row["Start_Position"]), str(row["Tumor_Seq_Allele2"])) in hotspot_set,
                axis=1
            )
        else:
            out["hotspot"] = False
        
        return out

    def mark_phased_vars(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate phasing:
          - Parse LPS tokens (comma-delimited).
          - Count occurrences; tokens with count >=2 define multi-member phase groups.
          - is_phased: True if variant belongs to any multi-member token.
          - phase_set_size: max count among its multi-member tokens (0 if none).
          - variant_source initialized / augmented with 'hotspot' and/or 'phase_group'.
        """
        out = df.copy()
        # Normalize/initialize variant_source as list-of-tags per row
        if "variant_source" not in out.columns:
            out["variant_source"] = [[] for _ in range(out.shape[0])]
        else:
            out["variant_source"] = out["variant_source"].apply(
                lambda v: v.split(",") if isinstance(v, str) and v else (v if isinstance(v, list) else [])
            )

        # Parse tokens from guaranteed LPS column
        lps_tokens = out["LPS"].fillna("").astype(str).str.split(",")
        flat = lps_tokens.explode().str.strip()
        valid = flat[
            (~flat.str.lower().isin(BAD_TOKENS)) &
            (flat != "")
        ]
        counts = valid.value_counts()
        multi_member = set(counts[counts >= 2].index)

        phase_set_size: List[int] = []
        is_phased_flags: List[bool] = []

        for i, toks in enumerate(lps_tokens):
            toks_clean = [t.strip() for t in toks if t and t.lower() not in BAD_TOKENS]
            member_tokens = [t for t in toks_clean if t in multi_member]

            if member_tokens:
                size = max(counts[t] for t in member_tokens)
                phase_set_size.append(int(size))
                is_phased_flags.append(True)
                if "phase_group" not in out.at[i, "variant_source"]:
                    out.at[i, "variant_source"].append("phase_group")
            else:
                phase_set_size.append(0)
                is_phased_flags.append(False)

            # Append hotspot tag if applicable
            if bool(out.at[i, "hotspot"]) and "hotspot" not in out.at[i, "variant_source"]:
                out.at[i, "variant_source"].append("hotspot")

        out["phase_set_size"] = phase_set_size
        out["is_phased"] = is_phased_flags
        return out

    def filter_gnomad(self, df: pd.DataFrame, threshold: float) -> pd.DataFrame:
        """Keep rows where all gnomAD AF pop columns < threshold (or NaN), OR hotspot OR is_phased."""
        out = df.copy()
        cols = [c for c in GNOMAD_COLS if c in out.columns]
        if not cols:
            return out
        gnomad_ok = (out[cols].lt(threshold) | out[cols].isna()).all(axis=1)
        exempt = out["hotspot"]
        keep = gnomad_ok | exempt
        return out.loc[keep].copy()

    def remove_blacklisted(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drop variants flagged as blacklisted (blacklist always overrides exemptions)."""
        if "blacklist" not in df.columns:
            return df.copy()
        return df.loc[df["blacklist"] == False].copy()

    def filter_vaf(self, df: pd.DataFrame, min_vaf: float) -> pd.DataFrame:
        """
        Append 'high_vaf' to variant_source for any variant with AF >= min_vaf.
        Remove only variants that are:
          - NOT hotspot
          - NOT phased
          - AF < min_vaf
        Hotspot/phased variants are retained regardless of AF.
        """
        out = df.copy()
        out["variant_source"] = out["variant_source"].apply(
            lambda v: v if isinstance(v, list) else (v.split(",") if isinstance(v, str) and v else [])
        )
        keep_mask = []
        for i, row in out.iterrows():
            af = float(row.get("AF", 0.0))
            hotspot = bool(row.get("hotspot", False))
            phased = bool(row.get("is_phased", False))
            sources = row["variant_source"]
            passes_vaf = af >= min_vaf
            if passes_vaf and "high_vaf" not in sources:
                sources.append("high_vaf")

            if hotspot or phased or passes_vaf:
                keep_mask.append(True)
            else:
                keep_mask.append(False)

            out.at[i, "variant_source"] = sources

        out = out.loc[keep_mask].copy()
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
        """
        Apply depth/alt/UMI/low-normal AF filters to non-exempt variants only.
        Exempt variants (hotspot OR is_phased) always pass.
        """
        required = {"t_depth","n_depth","t_alt_count","UMI_3_count","AF","hotspot","is_phased"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {sorted(missing)}")

        out = df.copy()
        exempt = out["hotspot"] | out["is_phased"]

        # Base thresholds (non-exempt rows)
        base_mask = (
            (out["t_alt_count"] >= min_alt_tum) &
            (out["n_depth"] >= min_germline_depth) &
            (out["t_depth"] >= min_t_depth)
        )
        # Low-normal safeguard
        base_mask &= ~((out["n_depth"] <= low_normal_depth) & (out["AF"] > low_normal_AF))
        # UMI tier
        umi_mask = (
            ((out["t_alt_count"] >= low_alt_thresh) & (out["UMI_3_count"] >= min_UMI_3_count)) |
            ((out["t_alt_count"] < low_alt_thresh) & (out["UMI_3_count"] >= low_alt_UMI_3_count))
        )
        base_mask &= umi_mask

        final_mask = base_mask | exempt
        filtered = out.loc[final_mask].copy()
        return filtered

    def filter_background_rates(
        self,
        df: pd.DataFrame,
        background_file: str,
        min_samples: int = 20,
        n_std: float = 2.0
    ) -> pd.DataFrame:
        """
        Filter variants based on background mutation rates.
        
        - If background_file is empty or doesn't exist, skip filter
        - SNVs: check alt allele at position
        - INS/DEL: check INS/DEL rate at position
        - DNP/TNP/ONP: check all affected positions pass filter
        - Variants not in background index are removed (suspicious)
        - Positions with n_samples < min_samples are removed (insufficient data)
        - VAF must be >= (mean + n_std*std) for the change
        - Hotspot OR phased variants bypass all checks
        """
        # Skip if no background file provided
        if not background_file or not os.path.isfile(background_file):
            print("No background rates file provided or file not found, skipping background filter")
            return df.copy()
        
        out = df.copy()
        
        # Load background rates
        print(f"Loading background mutation rates from {background_file}...")
        bg = pd.read_csv(background_file, sep="\t")
        bg = bg.set_index(['chromosome', 'position'])
        print(f"Loaded {len(bg)} positions from background index")
        
        # Separate exempt variants (hotspot OR phased)
        exempt_mask = out["hotspot"] | out["is_phased"]
        exempt_df = out.loc[exempt_mask].copy()
        non_exempt_df = out.loc[~exempt_mask].copy()
        
        print(f"  - {len(exempt_df)} exempt variants (hotspot or phased)")
        print(f"  - {len(non_exempt_df)} non-exempt variants to filter")
        
        # Mark noisy variants in non-exempt set
        non_exempt_df["noisy"] = False
        
        for variant in non_exempt_df.itertuples():
            chrom = str(variant.Chromosome)
            pos = int(variant.Start_Position)
            ref_allele = str(variant.Reference_Allele)
            alt_allele = str(variant.Tumor_Seq_Allele2)
            vaf = float(variant.AF)
            var_type = str(variant.Variant_Type)
            
            # Determine positions and alt bases to check
            if var_type in ["DNP", "TNP", "ONP"]:
                # Multi-nucleotide: check all affected positions
                positions_to_check = list(range(pos, pos + len(ref_allele)))
                alt_bases = list(alt_allele)
            else:
                # SNV, INS, DEL: single position
                positions_to_check = [pos]
                alt_bases = [alt_allele]
            
            # Check each position - any failure marks variant as noisy
            for check_pos, check_alt in zip(positions_to_check, alt_bases):
                bg_key = (chrom, check_pos)
                
                # Position not in background
                if bg_key not in bg.index:
                    non_exempt_df.at[variant.Index, "noisy"] = True
                    break
                
                bg_row = bg.loc[bg_key]
                
                # Insufficient samples
                if bg_row['n_samples'] < min_samples:
                    non_exempt_df.at[variant.Index, "noisy"] = True
                    break
                
                # Determine column based on variant type
                if var_type == "INS":
                    mean_col, std_col = "INS_mean", "INS_std"
                elif var_type == "DEL":
                    mean_col, std_col = "DEL_mean", "DEL_std"
                else:
                    # SNV or MNP: use alt base
                    mean_col, std_col = f"{check_alt}_mean", f"{check_alt}_std"
                
                # Column missing - can't filter, assume pass
                if mean_col not in bg_row.index or std_col not in bg_row.index:
                    continue
                
                # Calculate threshold and check
                bg_mean = float(bg_row[mean_col])
                bg_std = float(bg_row[std_col])
                threshold = bg_mean + (n_std * bg_std)
                
                if vaf < threshold:
                    non_exempt_df.at[variant.Index, "noisy"] = True
                    break
        
        # Filter out noisy variants
        passed_df = non_exempt_df.loc[~non_exempt_df["noisy"]].drop(columns=["noisy"])
        
        # Report stats
        n_noisy = non_exempt_df["noisy"].sum()
        print(f"Background rate filtering: {n_noisy} variants marked as noisy and removed")
        
        # Combine passed variants with exempt variants
        result = pd.concat([passed_df, exempt_df], ignore_index=True)
        print(f"  - Total kept: {len(result)}/{len(out)} ({len(passed_df)} passed + {len(exempt_df)} exempt)")
        
        return result


def main():
    """Parse args, run the pipeline, and write the filtered MAF."""
    args = get_args()
    cfg = PipelineConfig(**vars(args))
    pipeline = VariantFilterPipeline(cfg)
    df = pipeline.filter_variants()
    print(f"Finished custom filtering, {df.shape[0]} variants remain")
    df.to_csv(cfg.output_maf, sep="\t", index=False)

if __name__ == '__main__':
    main()
