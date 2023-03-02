#!/usr/bin/env python3

import os
import warnings
import argparse
import numpy as np
import pandas as pd
import oncopipe as op
import math

def main():

    input_bam = snakemake.input[0]
    input_bai = snakemake.input[1]
    input_maf = open(snakemake.input[2], "r")
    outfile = open(snakemake.output[0], "w")

    # Skip sample if no variants in filtered MAF file
    line_count = 0
    for line in input_maf:
        line_count += 1
        if line_count > 1:
            break
    if line_count < 2:
        line = generate_igv_batch_footer()
        output_lines(line, outfile)
        input_maf.close()
        outfile.close()
        exit()

    # Return to top of MAF
    input_maf.seek(0)

    # Read MAF file containing variants and create a dataframe linking regions, sample_ids, and bam paths
    regions = get_regions_df(
        input_maf,
        seq_type=snakemake.params[2],
        padding=snakemake.params[3]
    )

    input_maf.close()

    # Format and output the batch script 
    generate_igv_batch(
        bam = input_bam,
        bai = input_bai,
        regions = regions,
        output = outfile,
        max_height = snakemake.params[4],
        seq_type = snakemake.params[2],
        genome_build = snakemake.params[1],
        snapshot_dir = snakemake.params[0],
        igv_options = snakemake.params[5],
        image_format = snakemake.params[6]
    )

    outfile.close()

def get_regions_df(input_maf, seq_type, padding):
    # Read MAF as dataframe
    maf = pd.read_table(input_maf, comment="#", sep="\t")

    # Make sure required minimum columns are present in the maf
    columns = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Hugo_Symbol",
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

    # Specify the regions that will be captured by IGV based on variant positions
    region_position = (maf["Start_Position"]).astype(str)
    snapshot_start = (maf["Start_Position"] - padding).astype(str)
    snapshot_end = (maf["End_Position"] + padding).astype(str)
    snapshot_coordinates = "chr" + chrom + ":" + snapshot_start + "-" + snapshot_end
    regions = "chr" + chrom + ":" + region_position

    regions_df = pd.DataFrame(
        {"chromosome": "chr" + chrom,
        "region": regions,
        "region_name": maf.Hugo_Symbol,
        "sample_id": maf.Tumor_Sample_Barcode,
        "snapshot_coordinates": snapshot_coordinates,
        "padding": padding
        }
    )
    
    return regions_df

def generate_igv_batch_header(bam_file, index_file, max_height, genome_build):
    lines = []

    genome_build = genome_build.replace("grch37","hg19").replace("grch38","hg38")

    bam_file = os.path.realpath(bam_file)
    lines.append(f"load {bam_file}")

    bai_file = os.path.realpath(index_file)
    lines.append(f"index={bai_file}")

    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"genome {genome_build}")

    return lines

def generate_igv_batch_per_row(coordinates, snapshot_filename, igv_options):
    lines = []
    lines.append(f"goto {coordinates}")
    lines.append("sort")
    lines.append("collapse")
    for option in igv_options:
        lines.append(option)
    lines.append(f"snapshot {snapshot_filename}")

    return lines

def generate_igv_batch(bam, bai, regions, output, max_height, seq_type, genome_build, snapshot_dir, igv_options, image_format):

    # Lines for batch script encompassing all regions and sample_ids
    all_lines = []

    header = generate_igv_batch_header(
        bam, bai, max_height, genome_build
    )

    all_lines.extend(header)

    for chrom in regions.chromosome.unique():
        chrom_regions = regions[regions["chromosome"]==chrom]

        lines = generate_igv_batch_per_region(
            regions=chrom_regions,
            max_height=max_height,
            seq_type=seq_type,
            genome_build=genome_build,
            snapshot_dir=snapshot_dir,
            options=igv_options,
            image_format=image_format
        )

        if lines is not None:
            all_lines.extend(lines)

    footer = generate_igv_batch_footer()
    all_lines.extend(footer)

    output_lines(all_lines, output)


def generate_igv_batch_per_region(regions, max_height, seq_type, genome_build, snapshot_dir, options, image_format):
    
    # Batch script lines
    lines = []

    # Set up snapshot directory string
    dir_chrom = regions.chromosome.unique()[0].split(":")[0]
    seq_type_build = f"{seq_type}--{genome_build}"

    # Add snapshot directory line to batch script
    snapshot_regions_dir = os.path.join(snapshot_dir, seq_type_build, dir_chrom, "")
    lines.append(f"snapshotDirectory {snapshot_regions_dir}")

    # Add lines to batch script for each sample
    for _, row in regions.iterrows():
        # Add components of filename as a list
        filename = []

        filename.append(row.region)

        filename.append(str(row.padding))

        filename.append(row.region_name)
        
        filename.append(row.sample_id)

        #if not image_format.startswith("."):
        #    image_format = "." + image_format

        filename = "--".join(filename) + ".png"

        row_lines = generate_igv_batch_per_row(coordinates = row.snapshot_coordinates, snapshot_filename = filename, igv_options = options)

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

if __name__ == "__main__":
    main()

