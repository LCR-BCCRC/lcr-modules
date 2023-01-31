#!/usr/bin/env python3

import os
import warnings
import argparse
import numpy as np
import pandas as pd
import oncopipe as op
import math

def main():
    # Parse arguments
    args = parse_arguments()

    # Read MAF file containing variants and create a dataframe linking regions, sample_ids, and bam paths
    regions = get_regions_df(
        args.input_maf,
        metadata=args.metadata,
        seq_type=args.seq_type,
        padding=args.padding)

    # Format and output the batch script 
    generate_igv_batch(
        regions = regions,
        output = args.output,
        max_height = args.max_height,
        seq_type = args.seq_type,
        genome_build = args.genome_build,
        snapshot_dir=args.snapshot_dir,
        n_snapshots=args.n_snapshots)

    close_files(args)

def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "input_maf",
        type=argparse.FileType("r"),
        default="-",
        help=f"Input MAF. Can be '-' for stdin"
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="OUTPUT_FILE",
        type=argparse.FileType("w"),
        help="Output IGV batch script."
    )

    parser.add_argument(
        "--metadata",
        "-v",
        metavar="METADATA",
        type=argparse.FileType("r"),
        help="Metadata mapping sample IDs to BAM paths"
    )

    default_padding = 300
    parser.add_argument(
        "--padding",
        "-p",
        type=int,
        default=default_padding,
        help=(
            "Amount of padding added before and after each locus. "
            f"Default padding is {default_padding}"
        ),
    )

    default_max_height = 400
    parser.add_argument(
        "--max_height",
        "-m",
        type=int,
        default=default_max_height,
        help="Maximum panel height in IGV. Default max height is {default_max_height}"
    )

    parser.add_argument(
        "--snapshot_dir",
        "-d",
        required=True,
        help=(
            "Parent directory where {chromosome}/{region} subdirectories will be "
            "populated with IGV snapshots."
        )
    )

    parser.add_argument(
        "--n_snapshots",
        "-n",
        type=int,
        default=20,
        help=(
            "Maximum number of different snapshots for each position."
        )
    )

    parser.add_argument(
        "--genome_build",
        "-g",
        required=True,
        help="Specify IGV genome build for snapshots."
    )

    parser.add_argument(
        "--seq_type",
        "-s",
        required=True,
        type=str,
        help="Specify sequencing type for BAM extraction."
    )

    args = parser.parse_args()

    return args

def get_regions_df(input_maf, metadata, seq_type, padding):
    # Read MAF as dataframe
    maf = pd.read_table(input_maf, comment="#")

    # Read metadata as dataframe
    metadata = pd.read_table(metadata, comment="#")

    # Filter metadata down to only samples of required seq_type
    metadata = metadata[metadata["seq_type"]==seq_type]
    metadata = metadata[["sample_id","link_name"]]

    # Make sure required minimum columns are present in the maf
    columns = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Tumor_Sample_Barcode",
    ]
    
    assert(all(c in list(maf.columns) for c in columns)), (
        "The following required columns are missing: "
        f"{[columns[missing_ix] for missing_ix in [ix for ix, bool_val in enumerate([col not in list(maf.columns) for col in columns]) if bool_val==True]]}"
    )

    # Check if there are issues with input file
    for column in columns:
        column_values = maf[column]
        is_any_na = pd.isna(column_values).any()
        assert not is_any_na, (
            f"The '{column}' column contains NA values. This might be caused "
            "by an incorrectly formatted input MAF file. Please ensure that "
            f"all of the following columns have values: {', '.join(columns)}."
            f"Here's a preview of the MAF file after being parsed:\n\n {maf}"
        )

    # Create a pandas dataframe with to link regions with sample_ids and bam files

    chrom = (maf["Chromosome"].astype(str)).apply(lambda x: x.replace("chr",""))

    # Snapshots will be held in parent directories of 1000-nt intervals for easier navigation
    dir_start = ((maf["Start_Position"] / 1000).apply(lambda x: math.trunc(x)) * 1000).astype(str)
    dir_end = (dir_start.astype(int) + 1000).astype(str)
    dir_regions = "chr" + chrom + ":" + dir_start + "_" + dir_end

    # Specify the regions that will be captured by IGV based on variant positions and padding
    region_start = (maf["Start_Position"] - padding).astype(str)
    region_end = (maf["End_Position"] + padding).astype(str)
    regions = "chr" + chrom + ":" + region_start + "-" + region_end

    regions_df = pd.DataFrame(
        {"dir_regions": dir_regions,
        "regions": regions,
        "region_name": maf.Hugo_Symbol,
        "sample_id": maf.Tumor_Sample_Barcode,
        }
    )

    # Link bam paths to regions by merging metadata and regions dataframes by sample_id 
    regions_df = pd.merge(regions_df, metadata, on="sample_id", how="left")
    
    return regions_df

def generate_igv_batch_header(bam_file, max_height, snapshot_dir, genome_build):
    lines = []

    bam_file = os.path.realpath(bam_file)
    lines.append(f"load {bam_file}")

    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"snapshotDirectory {snapshot_dir}")
    lines.append(f"genome {genome_build}")

    return lines

def generate_igv_batch_per_row(regions, snapshot_filename):
    lines = []
    lines.append(f"goto {regions}")
    lines.append("sort")
    lines.append("collapse")
    lines.append(f"snapshot {snapshot_filename}")
    lines.append("new")

    return lines

def generate_igv_batch_per_region(regions, max_height, genome_build, snapshot_dir):
    
    # Lines of batch script
    lines = []

    # Add lines to batch script for each region
    for _, row in regions.iterrows():
        filename = []

        filename.append(row.regions)
        if "region_name" in row:
            filename.append(row.region_name)
        filename.append(row.sample_id)

        filename = "--".join(filename) + ".png"
        filename = filename.replace(" ", "_")

        bam_file = row.link_name

        genome_build = genome_build.replace("grch37","hg19").replace("grch38","hg38")

        header = generate_igv_batch_header(
            bam_file, max_height, snapshot_dir, genome_build
        )
        lines.extend(header)

        row_lines = generate_igv_batch_per_row(regions = row.regions, snapshot_filename = filename)

        lines.extend(row_lines)
    
    return lines

def close_files(args):
    args_dict = vars(args)
    for arg_value in args_dict.values():
        if hasattr(arg_value, "close"):
            arg_value.close()

def generate_igv_batch_footer():
    lines = []
    lines.append("exit")
    return lines

def output_lines(lines, output):
    lines.append("")
    text = "\n".join(lines)
    output.write(text)

def generate_igv_batch(regions, output, max_height, seq_type, genome_build, snapshot_dir, n_snapshots):

    # The lines for the batch script encompassing all regions and sample_ids
    all_lines = []

    # Create batch scripts per unique 1000-nt interval region
    dir_regions = regions.dir_regions.unique()

    for dir_region in dir_regions:

        # Get chromosome string for parent directory
        dir_chrom = dir_region.split(":")[0]
        # Get 1000nt interval region for subdirectory
        dir_interval = dir_region.split(":")[1]

        seq_type_build = f"{seq_type}--{genome_build}"
        
        region_snapshot_dir = os.path.join(snapshot_dir, seq_type_build, dir_chrom, dir_interval, "")
        
        # Subset all regions down to those in the 1000nt interval
        regions_in_dir = regions[regions["dir_regions"]==dir_region]

        # Iterate through unique regions within interval
        for unique_region in regions_in_dir.regions.unique():

            # Subset rows by number of snapshots desired per region
            regions_subset = regions_in_dir[regions_in_dir["regions"]==unique_region][:n_snapshots]
            
            # Generate lines of batch script per region subset
            lines = generate_igv_batch_per_region(
                regions = regions_subset,
                max_height=max_height, 
                genome_build=genome_build,
                snapshot_dir=region_snapshot_dir)

            if lines is not None:
                all_lines.extend(lines)

    footer = generate_igv_batch_footer()
    all_lines.extend(footer)

    output_lines(all_lines, output)

if __name__ == "__main__":
    main()

