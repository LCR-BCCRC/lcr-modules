"""Variants are filtered through multiple steps:

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
9) Minimum UMI support filter is applied (UMI_max => 3).
10) A variant_source column is added indicating: phase_group, high_vaf, hotspot, 
    or combinations thereof.
"""
import pandas as pd
from collections import Counter
import argparse


GNOMAD_COLS = ["gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"]


def get_args():
    parser = argparse.ArgumentParser()
    # io
    parser.add_argument('--input_maf',required=True,type=str,help='')
    parser.add_argument('--output_maf',required=True,type=str,help='')
    parser.add_argument('--blacklist',required=True,type=str,help='')
    parser.add_argument('--hotspots',required=True,type=str,help='')
    # filtering params
    parser.add_argument('--gnomad_threshold',required=True,type=float,help='')
    # hard filters
    parser.add_argument('--min_alt_depth_tum',required=True,type=int,help='')
    parser.add_argument('--min_germline_depth',required=True,type=int,help='')
    parser.add_argument('--min_tumour_vaf',required=True,type=float,default=0.01,help='')
    parser.add_argument('--min_t_depth',required=False,type=int,default=200,help='')
    # umi support filters
    parser.add_argument('--min_UMI_3_count',required=False,type=int,default=1,help='Minimum UMI_3_count for variants to pass filter,')
    parser.add_argument('--low_alt_thresh',required=False,type=int,default=20,help='If alt read count is below this threshold then min_alt_UMI_3_count is used instead of min_UMI_3_count')
    parser.add_argument('--low_alt_min_UMI_3_count',required=False,type=int,default=2,help='Minimum UMI_3_count for variants with low alt read count to pass filter. Allows requiring more UMI support.')
    # low normal depth filter, in cases where germlines vars may have been missed
    parser.add_argument('--low_normal_depth',required=False,type=int,default=100,help='If normal depth is below this threshold then variants with AF > low_normal_AF are removed')
    parser.add_argument('--low_normal_AF',required=False,type=float,default=0.3,help='If normal depth is below low_normal_depth then variants with AF > this threshold are removed')

    return parser.parse_args()

def read_maf(maf_file: str) -> pd.DataFrame:
    """Read a maf file into a pandas DataFrame.

    And add variant_key column that is: chr:start

    Args:
        maf_file (str): The path to the maf file.
    Returns:
        pd.DataFrame: The DataFrame of the maf file.
    """

    indf= pd.read_csv(maf_file, sep="\t")
    indf["variant_key"] = indf["Chromosome"].astype(str) + ":" + indf["Start_Position"].astype(str)
    return indf.copy()

def mark_blacklist_hotspot(df: pd.DataFrame, blacklist: str, hotspots: str) -> pd.DataFrame:
    """Mark variants that are in the blacklist or hotspots file.

    Args:
        df (pd.DataFrame): DataFrame of variants.
        blacklist (str): Path to the blacklist file.
        hotspots (str): Path to the hotspots file.
    Returns:
        pd.DataFrame: DataFrame of variants with a new column "blacklist_hotspot" that is True if the variant is in the blacklist or hotspots file.
    """
    # read in blacklist and hotspots
    blacklist_df = pd.read_csv(blacklist, sep="\t")
    blacklist_df.columns = ["variant_key"]
    hotspots_df = pd.read_csv(hotspots, sep="\t")
    hotspots_df.columns = ["variant_key"]

    # if variant_key in black list mark hotspot as True
    df["blacklist"] = df["variant_key"].isin(blacklist_df["variant_key"])
    df["hotspot"] = df["variant_key"].isin(hotspots_df["variant_key"])

    return df.copy()

def filter_gnomad(df : pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Filter variants that have a gnomadAF below a certain threshold.

        Ignore NAs in the gnomadAF columns.
    Args:
        df (pd.DataFrame): DataFrame of variants.
        threshold (float): The threshold to filter variants by.
    Returns:
        pd.DataFrame: DataFrame of variants that have a gnomadAF below the threshold.
    """
    for col in GNOMAD_COLS:
        # is less than or is na
        df = df[(df[col] < threshold) | (df[col].isna())]
    return df.copy()

def recalculate_af(df:pd.DataFrame) -> pd.DataFrame:
    """Recalculate the AF column based on the new set of variants.

    Args:
        df (pd.DataFrame): DataFrame of variants.
    Returns:
        pd.DataFrame: DataFrame of variants with recalculated AF column.
    """
    # calculate the AF
    df["AF"] = df["t_alt_count"] / df["t_depth"]
    return df.copy()

def _fetch_phase_sets(indf: pd.DataFrame, bad_tokens: list) -> set:
    """Fetch phase set memberships for each variant from LPS column
    """

    lps_col = indf["LPS"].fillna("").astype(str).str.split(",").explode().str.strip()
    
    phase_sets = lps_col[~lps_col.str.lower().isin(bad_tokens)]
    # summary count phase sets
    phase_set_counts = Counter(phase_sets)
    
    # find phase sets with more than 1 variant
    multi_variant_phase_sets = {ps for ps, count in phase_set_counts.items() if count > 1}

    if len(multi_variant_phase_sets) == 0:
        multi_variant_phase_sets = None

    return multi_variant_phase_sets


def filter_vaf_and_phase(indf: pd.DataFrame, min_vaf: float) -> pd.DataFrame:
    """Filter variants that are below min VAF unless they
    are a part of a phase set of variants that also made it
    this far in the filtering process.

    Args:
        indf (pd.DataFrame): DataFrame of variants.
    Returns:
        pd.DataFrame: DataFrame of variants that pass filters.
    """
    BAD_TOKENS = {"", "na", "n/a", "none", "nan"}
    outdf = indf.copy()
    outdf["variant_source"] = [[] for _ in range(outdf.shape[0])]
    # get count of phase sets in LPS col
    multi_var_phase_sets = _fetch_phase_sets(outdf, BAD_TOKENS)
    print(f"Phase sets to keep: {multi_var_phase_sets}")

    for row in outdf.itertuples():
        lsps = [t.strip() for t in str(row.LPS).split(",")]
        lsps = [t for t in lsps if t and t.lower() not in BAD_TOKENS]

        if multi_var_phase_sets and any(ps in multi_var_phase_sets for ps in lsps):
            row.variant_source.append("phase_group")
        if row.AF >= min_vaf:
            row.variant_source.append("high_vaf")
        if bool(row.hotspot):
            row.variant_source.append("hotspot")

    # only keep variants with a len(variant_source) > 0
    outdf["variant_source"] = outdf["variant_source"].apply(lambda x: ",".join(x) if len(x) > 0 else None)
    # remove variants that didnt get any filter label
    outdf = outdf[outdf["variant_source"].notna()].copy()

    return outdf

def mark_potential_chip(indf: pd.DataFrame) -> pd.DataFrame:
    """Mark potential CHIP mutations.

    A two part filter:
    1) VAF in tumour has to be 3x that in normal

    2) Any variant with more than 5 reads alt support in normal is removed

    If a variant fails either of these filters it is marked as CHIP.
    """
    chip_genes = ["DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "JAK2", "SF3B1", "SRSF2"]

    indf = indf.copy()
    # calc VAF of normal
    indf["VAF_normal"] = indf["n_alt_count"] / indf["n_depth"]
    # calc var vaf (gets messed up in original maf due to augment_ssm)
    indf["AF"] = indf["t_alt_count"] / indf["t_depth"]

    for variant in indf.itertuples():
        if variant.Hugo_Symbol not in chip_genes:
            try:
                if (variant.VAF_normal / variant.AF) >= 3 or (variant.VAF_normal >= 0.02):
                    indf.at[variant.Index, "CHIP"] = True
                else:
                    indf.at[variant.Index, "CHIP"] = False
            except ZeroDivisionError:
                print(f"Error in calculating CHIP for {variant.Index} in {indf}")
        # be more strict with common CHIP genes, if 1 read in normal then mark as CHIP
        else:
            try:
                if (variant.VAF_normal / variant.AF) >= 3 or (variant.VAF_normal >= 0.02) or (variant.n_alt_count >= 1):
                    indf.at[variant.Index, "CHIP"] = True
                else:
                    indf.at[variant.Index, "CHIP"] = False
            except ZeroDivisionError:

                print(f"Error in calculating CHIP for {variant.Index} in {indf}")
    return indf.copy()

def filter_by_read_support(
        indf: pd.DataFrame,
        min_t_depth: int,
        min_germline_depth: int,
        min_alt_tum: int,
        min_UMI_3_count: int = 1,

        low_alt_thresh: int = 20,
        low_alt_UMI_3_count: int = 2,
        
        low_normal_depth: int = 100,
        low_normal_AF: float = 0.3
        ) -> pd.DataFrame:
    """Filter variants by read support.
    A three part filter:
    1) Variants must have at least min_alt_tum in the tumour and min_germline_depth in the germline.
    2) Variants with low normal depth (<=100) and high AF (>0.3) are removed.
    3) Variants with low alt read count (<20) must have at least min_alt_UMI_3_count UMI_3_count, otherwise
       they must have at least min_UMI_3_count UMI_3_count.

    Hotspots are kept regardless of read support.

    """
    # check for required cols
    cols = set(indf.columns)
    req_cols = {"t_depth", "n_depth", "t_alt_count", "UMI_3_count", "AF"}
    if not req_cols.issubset(cols):
        missing = req_cols - cols
        raise ValueError(f"Input DataFrame is missing required columns: {missing}")
    
    # hard filters, create a mask series for each variant as boolean
    mask = (indf["t_alt_count"] >= min_alt_tum) & (indf["n_depth"] >= min_germline_depth) & (indf["t_depth"] >= min_t_depth)

    # but in cases with low n_depth be stricter on AF, in case any germlines werent properly filtered
    # due to low cov in the matched normal
    mask &= ~((indf["n_depth"] <= low_normal_depth) & (indf["AF"] > low_normal_AF))

    # apply UMI filtering, based on alt read count threshold
    umi_mask = (
                (indf["t_alt_count"] >= low_alt_thresh) & (indf["UMI_3_count"] >= min_UMI_3_count) |
                (indf["t_alt_count"] < low_alt_thresh) & (indf["UMI_3_count"] >= low_alt_UMI_3_count)
                )

    mask &= umi_mask

    # make sure all hotspots are kept
    mask |= (indf["hotspot"] == True)

    return indf[mask].copy()

def main():
    args = get_args()
    print("Starting custom filtering")

    # read input maf
    inmaf = read_maf(args.input_maf)
    print('\033[94m' + f"Input maf: {args.input_maf} with {inmaf.shape[0]} variants")
    # recalculate AF
    inmaf = recalculate_af(inmaf)
    # filter gnomad
    inmaf = filter_gnomad(inmaf, args.gnomad_threshold)
    # mark chip mutations
    inmaf = mark_potential_chip(inmaf)
    # mark blacklist and hotspot
    inmaf = mark_blacklist_hotspot(inmaf, args.blacklist, args.hotspots)
    # remove all blacklisted vars
    inmaf = inmaf.loc[inmaf["blacklist"]== False].copy()
    # filter vaf and phase
    outmaf = filter_vaf_and_phase(inmaf, args.min_tumour_vaf)

    # filter by read support
    inmaf = filter_by_read_support(
        outmaf,
        min_t_depth=args.min_t_depth,
        min_germline_depth=args.min_germline_depth,
        min_alt_tum=args.min_alt_depth_tum,
        min_UMI_3_count=args.min_UMI_3_count,
        low_alt_thresh=args.low_alt_thresh,
        low_alt_UMI_3_count=args.low_alt_min_UMI_3_count,
        low_normal_depth=args.low_normal_depth,
        low_normal_AF=args.low_normal_AF
    )

    print(f"Finished custom filtering, {outmaf.shape[0]} variants remain")
    # write output maf
    outmaf.to_csv(args.output_maf, sep="\t", index=False)

if __name__ == '__main__':
    main()