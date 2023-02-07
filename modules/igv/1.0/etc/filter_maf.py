#!/usr/bin/env python

import os
import math
import pandas as pd
import oncopipe as op

def filter_by_bed(maf, regions, metadata):

    # Remove rows that contain column names
    regions = regions[regions[0].str.contains("chrom")==False]

    # Create common columns between BED and MAF
    regions["chr_std"] = regions.apply(lambda x: str(x[0]).replace("chr",""), axis=1)
    regions["genomic_pos_std"] = regions["chr_std"] + ":" + regions[1].map(str)

    maf["chr_std"] = maf.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
    maf["genomic_pos_std"] = maf["chr_std"] + ":" + maf["Start_Position"].map(str)

    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def filter_by_maf(maf, regions, metadata):

    # Create common column by which to subset MAF
    for df in [maf, regions]:
        df["chr_std"] = df.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
        df["genomic_pos_std"] = df["chr_std"] + ":" + df["Start_Position"].map(str) + "_" + df["End_Position"].map(str)

    # Subset the MAF
    filtered_maf = maf[maf["genomic_pos_std"].isin(regions["genomic_pos_std"])]
    return filtered_maf

def maf_filter(maf, regions, regions_format, metadata, genome_build, seq_type, genome_map):
    # Read input MAF and regions file as dataframes
    maf_df = pd.read_table(maf, comment="#", sep="\t")

    if regions_format != "bed":
        regions_df = pd.read_table(regions, comment="#", sep="\t")
    else:
        regions_df = pd.read_table(regions, comment="#", sep="\t", header=None)

    # Select rows in MAF containing correct seq_type and build
    metadata = op.filter_samples(metadata, seq_type=seq_type)
    genome_build_list = genome_map[genome_build]
    metadata = op.filter_samples(metadata, genome_build=genome_build_list)

    maf_df = maf_df[maf_df["Tumor_Sample_Barcode"].isin(metadata.sample_id)]

    filter_functions = {
        "maf": filter_by_maf,
        "bed": filter_by_bed
        }
    
    return filter_functions[regions_format](maf_df, regions_df, metadata)

def maf_reduce_snapshots(maf, snapshots):
    # Only include max of number of snapshots for each variant
    maf = maf.groupby(["Chromosome","Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]).head(n=snapshots)

    return maf

def write_output(maf, outfile):
    maf.to_csv(outfile, sep="\t", index=False)

maf_file = snakemake.input[0]

regions_file = snakemake.input[1]
regions_format = snakemake.params[0]

if regions_format == "oncodriveclustl":
    # This should act as a global variable
    CLUSTL_PARAMS = snakemake.params[5]

n_snapshots = snakemake.params[6]

# Metadata file or dataframe mapping sample_ids to bam file paths
metadata = snakemake.params[1]
if not isinstance(metadata, pd.DataFrame):
    metadata = pd.read_table(metadata, comment="#", sep="\t")

maf_genome_build = snakemake.params[2]
maf_seq_type = snakemake.params[3]

# Dictionary of genome builds present in the MAF to label as grch37 / hg38
genome_map = snakemake.params[4]

output_file = snakemake.output[0]

# Peform filtering

filtered_maf = maf_filter(
    maf=maf_file, 
    regions=regions_file,
    regions_format=regions_format,
    metadata=metadata,
    genome_build=maf_genome_build,
    seq_type=maf_seq_type,
    genome_map=genome_map
    )

filtered_maf = maf_reduce_snapshots(maf=filtered_maf, snapshots=n_snapshots)

write_output(filtered_maf, output_file)
