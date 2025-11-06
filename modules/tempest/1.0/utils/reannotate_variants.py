"""
Script to reannoate variants with the target names in panel bed file.

The reference data used for annotation is not always the same that was
used to create the panel. 

Example usage:

python reannotate_variants.py --panel_bed panel.bed --in_maf variants.maf --out_maf variants_reannotated.maf

Arguments:
    --panel_bed: Path to the panel bed file.
    --in_maf: Path to the variants maf file.
    --out_maf: Path to the output maf file.

tier list must have the columns "Gene" and "Tier"

"""
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Reannotate variants with the target names in panel bed file.")
    parser.add_argument("--panel_bed", type=str, required=True, help="Path to the panel bed file.")
    parser.add_argument("--in_maf", type=str, required=True, help="Path to the variants maf file.")
    parser.add_argument("--out_maf", type=str, required=True, help="Path to the output maf file.")
    return parser.parse_args()

def read_panel_bed(panel_bed: str) -> pd.DataFrame:
    """
    Read the panel bed file and return a dictionary of target names and their coordinates.

    Args:
        panel_bed: Path to the panel bed file.

    Returns:
        A pandas dataframe with the target names and their coordinates.
    """
    panel_df = pd.read_csv(panel_bed, sep="\t", header=None)
    panel_df.columns = ["chrom", "start", "end", "name"]
    # make sure start and end are integers
    panel_df["start"] = panel_df["start"].astype(int)
    panel_df["end"] = panel_df["end"].astype(int)
    return panel_df

def read_maf(maf: str) -> pd.DataFrame:
    """
    Read the maf file and return a pandas dataframe.

    Args:
        maf: Path to the maf file.

    Returns:
        A pandas dataframe with the maf file.
    """
    maf_df = pd.read_csv(maf, sep="\t")
    # make sure start and end are integers
    maf_df["Start_Position"] = maf_df["Start_Position"].astype(int)
    maf_df["End_Position"] = maf_df["End_Position"].astype(int)
    return maf_df

def reannotate_variants(maf_df: pd.DataFrame, panel_df: pd.DataFrame) -> pd.DataFrame:
    """
    Reannotate the variants with the target names in panel bed file.

    Search the bed file for the target name in the maf file, and change
    the Hugo_Symbol column to the target name.

    Args:
        maf_df: A pandas dataframe with the maf file.
        panel_df: A pandas dataframe with the target names and their coordinates.

    Returns:
        A pandas dataframe with the reannotated maf file.
    """
    indf = maf_df.copy()
    for index, row in indf.iterrows():
        try:
            target_name = panel_df.loc[(panel_df["chrom"] == row["Chromosome"]) & (panel_df["start"] <= row["Start_Position"]) & (panel_df["end"] >= row["End_Position"])]["name"].values[0]
            indf.at[index, "Hugo_Symbol"] = target_name
        except IndexError:
            print(f"No target name found for {row['Hugo_Symbol']} at {row['Chromosome']}:{row['Start_Position']}-{row['End_Position']}")
            continue
    return indf

def main():
    args = parse_args()
    panel_df = read_panel_bed(args.panel_bed)
    maf_df = read_maf(args.in_maf)
    maf_df = reannotate_variants(maf_df, panel_df)
    maf_df.to_csv(args.out_maf, sep="\t", index=False)


if __name__ == "__main__":
    main()