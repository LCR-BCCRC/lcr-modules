#!/usr/bin/env python

import os
import math
import pandas as pd
import oncopipe as op

def filter_by_bed(maf, regions):

    # Remove row containing column names
    regions = regions[regions[0].str.contains("chrom")==False]

    # Create common columns between BED and MAF
    regions["chr_std"] = regions.apply(lambda x: str(x[0]).replace("chr",""), axis=1)
    regions["genomic_pos_std"] = regions["chr_std"] + ":" + regions[1].map(str)

    maf["chr_std"] = maf.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
    maf["genomic_pos_std"] = maf["chr_std"] + ":" + maf["Start_Position"].map(str)

    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def filter_by_maf(maf, regions):

    # Create common column by which to subset MAF
    for df in [maf, regions]:
        df["chr_std"] = df.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
        df["genomic_pos_std"] = df["chr_std"] + ":" + df["Start_Position"].map(str)

    # Subset the MAF
    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def maf_filter(maf, regions, regions_format):
    
    if regions_format != "bed":
        regions_df = pd.read_table(regions, comment="#", sep="\t")
    else:
        regions_df = pd.read_table(regions, comment="#", sep="\t", header=None)

    # Return empty dataframe without filtering if df is empty
    if len(maf)==0:
        return maf

    filter_functions = {
        "maf": filter_by_maf,
        "bed": filter_by_bed
        }
    
    return filter_functions[regions_format](maf, regions_df)

def maf_reduce_snapshots(maf, snapshots):
    # Only include max of number of snapshots for each variant
    maf = maf.groupby(["Chromosome","Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]).head(n=snapshots)

    return maf

def maf_add_columns(maf, metadata):
    # Read input MAF as df
    maf = pd.read_table(maf, comment="#", sep="\t")

    sample_id = maf["Tumor_Sample_Barcode"].unique()[0]

    row = metadata[metadata["tumour_sample_id"]==sample_id]

    seq_type = row["tumour_seq_type"].item()
    genome_build = row["tumour_genome_build"].item()
    normal_sample_id = row["normal_sample_id"].item()
    pair_status = row["pair_status"].item()

    maf["seq_type"] = seq_type
    maf["genome_build"] = genome_build
    maf["normal_sample_id"] = normal_sample_id
    maf["pair_status"] = pair_status
    
    return maf

def write_output(maf, outfile):
    maf.to_csv(outfile, sep="\t", index=False)

maf_file = snakemake.input[0]

regions_file = snakemake.input[1]
regions_format = snakemake.params[0]

metadata = snakemake.params[3]

if regions_format == "oncodriveclustl":
    # This should act as a global variable
    CLUSTL_PARAMS = snakemake.params[1]

n_snapshots = snakemake.params[2]

output_file = snakemake.output[0]

maf = maf_add_columns(maf=maf_file, metadata=metadata)

# Peform filtering

filtered_maf = maf_filter(
    maf=maf, 
    regions=regions_file,
    regions_format=regions_format
    )

filtered_maf = maf_reduce_snapshots(maf=filtered_maf, snapshots=n_snapshots)

write_output(filtered_maf, output_file)
