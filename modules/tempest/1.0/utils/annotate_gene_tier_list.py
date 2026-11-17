"""
Script to annotate the gene tier from the gene tier list
database provided.

Example tier list
/home/kyakimovich/repos/LLMPP/resources/curated/dlbcl_genes.tsv

"""
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate the gene tier from the gene tier list database provided.")
    parser.add_argument("--gene_tier_list", type=str, required=True, help="Path to the gene tier list database.")
    parser.add_argument("--in_maf", type=str, required=True, help="Path to the variants maf file.")
    parser.add_argument("--out_maf", type=str, required=True, help="Path to the output maf file.")
    return parser.parse_args()


def read_maf(maf: str) -> pd.DataFrame:
    """
    Read the maf file and return a pandas dataframe.
    """
    maf_df = pd.read_csv(maf, sep="\t")
    return maf_df

def read_gene_tier_list(gene_tier_list: str) -> pd.DataFrame:
    """
    Read the gene tier list database and return a pandas dataframe.
    """
    gene_tier_df = pd.read_csv(gene_tier_list, sep="\t")
    return gene_tier_df

def annotate_gene_tier(maf_df: pd.DataFrame, gene_tier_df: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate the gene tier from the gene tier list database provided.

    Need to accomodate not exact matches, so look to see if the gene from the tier list
    appears in the Hugo_Symbol column of the maf file.

    Args:
        maf_df: A pandas dataframe with the maf file.
        gene_tier_df: A pandas dataframe with the gene tier list database.

    Returns:
        A pandas dataframe with the annotated maf file.
    """
    # make gene tier dict
    gene_tier_dict = dict(zip(gene_tier_df["Gene"], gene_tier_df["Tier"]))

    for index, row in maf_df.iterrows():
        for gene in gene_tier_dict.keys():
            if gene in row["Hugo_Symbol"]:
                maf_df.at[index, "Tier"] = gene_tier_dict[gene]
                break
    return maf_df

def main():
    args = parse_args()
    maf_df = read_maf(args.in_maf)
    gene_tier_df = read_gene_tier_list(args.gene_tier_list)
    ann_maf_df = annotate_gene_tier(maf_df, gene_tier_df)
    ann_maf_df.to_csv(args.out_maf, sep="\t", index=False)

if __name__ == "__main__":
    main()
