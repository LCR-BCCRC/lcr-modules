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
import argparse

GNOMAD_COLS = ["gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"]


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_maf',required=True,type=str,help='')
    parser.add_argument('--output_maf',required=True,type=str,help='')
    parser.add_argument('--blacklist',required=True,type=str,help='')
    parser.add_argument('--hotspots',required=True,type=str,help='')
    parser.add_argument('--gnomad_threshold',required=True,type=float,help='')
    parser.add_argument('--min_alt_depth_tum',required=True,type=int,help='')
    parser.add_argument('--min_germline_depth',required=True,type=int,help='')
    parser.add_argument('--min_tumour_vaf',required=True,type=float,default=0.01,help='')
    parser.add_argument('--min_t_depth',required=False,type=int,default=200,help='')
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

def filter_vaf_and_phase(indf: pd.DataFrame, min_vaf: float) -> pd.DataFrame:
    """Filter variants that are below 0.1% VAF unless they
    are a part of a phase set of variants that also made it
    this far in the filtering process.

    Args:
        indf (pd.DataFrame): DataFrame of variants.
    Returns:
        pd.DataFrame: DataFrame of variants that are below 0.01 VAF (1%) unless they
        are a part of a phase set of variants that also made it
        this far in the filtering process. Or in a hotspot.
    """
    # get count of phase sets in LPS col
    phase_set_counts = indf["LPS"].value_counts().reset_index()
    # get phase sets that have more than 1 variant
    phase_sets_to_keep = phase_set_counts[(phase_set_counts["count"] > 1)].copy()
    print(f"Phase sets to keep: {phase_sets_to_keep}")

    outdf = indf.loc[indf["LPS"].isin(phase_sets_to_keep["LPS"]) | (indf["AF"] >= min_vaf) | (indf["hotspot"] == True)].reset_index(drop=True).copy()

    # label variant_source and create column
    if outdf.shape[0] > 0:
        for variant in outdf.itertuples():
            source_label= []
            if variant.LPS in phase_sets_to_keep["LPS"].tolist():
                source_label.append("phase_group")
            if variant.AF >= 0.01:
                source_label.append("high_vaf")
            if variant.hotspot == True:
                source_label.append("hotspot")
            if len(source_label) == 0:
                source_label.append("unknown")
            outdf.at[variant.Index, "variant_source"] = ",".join(source_label)
    else:
        outdf = outdf.assign(variant_source=None)

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
                if (variant.VAF_normal / variant.AF) > 3 or (variant.VAF_normal > 0.02):
                    indf.at[variant.Index, "CHIP"] = True
                else:
                    indf.at[variant.Index, "CHIP"] = False
            except ZeroDivisionError:
                print(f"Error in calculating CHIP for {variant.Index} in {indf}")
        # be more strict with common CHIP genes, if 1 read in normal then mark as CHIP
        else:
            try:
                if (variant.VAF_normal / variant.AF) > 3 or (variant.VAF_normal > 0.02) or (variant.n_alt_count > 1):
                    indf.at[variant.Index, "CHIP"] = True
                else:
                    indf.at[variant.Index, "CHIP"] = False
            except ZeroDivisionError:

                print(f"Error in calculating CHIP for {variant.Index} in {indf}")
    return indf.copy()


def min_read_support(indf: pd.DataFrame, min_alt_tum: int, min_germline_depth: int) -> pd.DataFrame:
    """
    Filter variants that have less than min_alt_tum in the 
    tumour or less than min_germline_depth in the germline.
    """
    outdf = indf[(indf["t_alt_count"] >= min_alt_tum) & (indf["n_depth"] >= min_germline_depth)].copy()

    # but in cases with low n_depth be stricter on AF, in case any germlines werent properly filtered
    # due to low cov in the matched normal
    outdf = outdf[~((outdf["n_depth"] <= 100) & (outdf["AF"] > 0.3))].copy()

    return outdf

def min_UMI_support(indf: pd.DataFrame, min_UMI: int,) -> pd.DataFrame:
    """Filter variants for a min UMI_3_count value, if the column exists
    
    or UMI_max
    It is a count of the numer of reads that have a UMI family of 3 or more.
    """

    if "UMI_max" in indf.columns:
        return indf[indf["UMI_max"] >= min_UMI].copy()
    else:
        return indf.copy()

def min_t_depth(indf: pd.DataFrame, min_t_depth: int) -> pd.DataFrame:
    """Filter variants for a min t_depth value, if the column exists

    Args:
        indf (pd.DataFrame): DataFrame of variants.
        min_t_depth (int): The minimum t_depth to filter variants by.
    Returns:
        pd.DataFrame: DataFrame of variants that have a t_depth greater than or equal to min_t_depth.
    """
    if "t_depth" in indf.columns:
        return indf[indf["t_depth"] >= min_t_depth].copy()
    else:
        return indf.copy()

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
    # filter chip mutations
    inmaf = mark_potential_chip(inmaf)
    # min alt support
    inmaf = min_read_support(inmaf, args.min_alt_depth_tum, args.min_germline_depth)
    # min t depth
    inmaf = min_t_depth(inmaf, args.min_t_depth)
    # mark blacklist and hotspot
    inmaf = mark_blacklist_hotspot(inmaf, args.blacklist, args.hotspots)
    # remove all blacklisted vars
    inmaf = inmaf.loc[inmaf["blacklist"]== False].copy()
    # filter vaf and phase
    outmaf = filter_vaf_and_phase(inmaf, args.min_tumour_vaf)
    # filter min UMI_max
    outmaf = min_UMI_support(outmaf, 3)
    print(f"Finished custom filtering, {outmaf.shape[0]} variants remain")
    # write output maf
    outmaf.to_csv(args.output_maf, sep="\t", index=False)

if __name__ == '__main__':
    main()