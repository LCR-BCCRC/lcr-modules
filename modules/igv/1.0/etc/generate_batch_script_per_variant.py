#!/usr/bin/env python3

import os
import warnings
import numpy as np
import pandas as pd
import oncopipe as op
import sys
import logging
import traceback

def log_exceptions(exctype, value, tb):
    logging.critical(''.join(traceback.format_tb(tb)))
    logging.critical('{0}: {1}'.format(exctype, value))

sys.excepthook = log_exceptions

def main():
    with open(snakemake.log[0], "w") as stdout:
        # Set up logging
        sys.stdout = stdout

        try:
            input_maf = open(snakemake.input[0], "r")
            input_bam = snakemake.input[1]
            input_bai = snakemake.input[2]

            # Skip if no variants in outfile
            line_count = 0
            for line in input_maf:
                line_count += 1
                if line_count > 1:
                    break
            if line_count < 2:
                input_maf.close()
                touch_output = open(snakemake.output[0],"w")
                touch_output.close()
                exit()

            # Return to top of MAF
            input_maf.seek(0)

            # Read MAF file and create dataframe
            regions = get_regions_df(
                input_maf,
                seq_type=snakemake.params[2],
                padding=snakemake.params[4]
            )

            input_maf.close()

            # Create the batch scripts
            generate_igv_batches(
                regions = regions,
                bam = input_bam,
                bai = input_bai,
                output_dir = snakemake.params[0],
                snapshot_dir = snakemake.params[1],
                genome_build = snakemake.params[2],
                seq_type = snakemake.params[3],
                igv_options = snakemake.params[5],
                max_height = snakemake.params[6],
                suffix = snakemake.params[7],
                as_pairs = snakemake.params[8],
                sleep_timer = snakemake.params[9]
            )

            touch_output = open(snakemake.output[0], "w")
            touch_output.close()
        
        except Exception as e:
            logging.error(e, exc_info=1)
            raise

def get_regions_df(input_maf, seq_type, padding):
    # Read MAF as dataframe
    maf = pd.read_table(input_maf, comment="#", sep="\t")

    chrom = (maf["Chromosome"].astype(str)).apply(lambda x: x.replace("chr",""))

    # Specify regions that will be captured by IGV based on variant positions
    region_start = (maf["Start_Position"]).astype(str)
    snapshot_start = (maf["Start_Position"] - padding).astype(str)
    snapshot_end = (maf["End_Position"] + padding).astype(str)
    snapshot_coordinates = "chr" + chrom + ":" + snapshot_start + "-" + snapshot_end
    regions = "chr" + chrom + ":" + region_start

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

def output_lines(lines, batch_output):
    output = open(batch_output, "w")
    lines.append("")
    text = "\n".join(lines)
    output.write(text)
    output.close()

def generate_igv_batch_per_row(coordinates, snapshot_filename):
    lines = []
    lines.append(f"goto {coordinates}")
    lines.append("sort")
    lines.append("collapse")
    lines.append(f"snapshot {snapshot_filename}")

    return lines

def generate_igv_batch_header(bam, index, max_height, genome_build, igv_options, sleep_timer, as_pairs):
    lines = []

    genome_build = genome_build.replace("grch37","hg19")

    bam_file = bam
    bai_file = index
    lines.append(f"load {bam_file} index={bai_file}")

    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"genome {genome_build}")
    
    if igv_options is not None:
        for option in igv_options:
            lines.append(option)

    lines.append(f"setSleepInterval {sleep_timer}")

    if as_pairs:
        lines.append("viewaspairs")

    return lines

def generate_igv_batches(regions, bam, bai, output_dir, snapshot_dir, genome_build, seq_type, igv_options, max_height, suffix, as_pairs=False, sleep_timer=2000):
    for _, row in regions.iterrows():
        all_lines = []

        header = generate_igv_batch_header(bam=bam, index=bai, max_height=max_height, genome_build=genome_build, igv_options=igv_options, sleep_timer=sleep_timer, as_pairs=as_pairs)
        all_lines.extend(header)

        dir_chrom = row.chromosome
        seq_type_build = f"{seq_type}--{genome_build}"

        snapshot_regions_dir = os.path.join(snapshot_dir, seq_type_build, dir_chrom, "")
        all_lines.append(f"snapshotDirectory {snapshot_regions_dir}")

        filename = []
        filename.append(row.region),
        filename.append(row.region_name)
        filename.append(row.sample_id)

        batch_filename = "--".join(filename) + suffix + ".batch"
        filename = "--".join(filename) + suffix + ".png"

        lines = generate_igv_batch_per_row(
            coordinates = row.snapshot_coordinates,
            snapshot_filename = filename
        )

        all_lines.extend(lines)

        for subdir in [os.path.join(output_dir, "single_batch_scripts"), os.path.join(output_dir, "single_batch_scripts", seq_type_build)]:
            if not os.path.exists(subdir):
                os.mkdir(subdir)

        batch_file_path = os.path.join(output_dir, "single_batch_scripts", seq_type_build, batch_filename)
        
        output_lines(all_lines, batch_file_path)

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[1],
        filemode='w'
    )
    main()
