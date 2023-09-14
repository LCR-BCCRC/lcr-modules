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
            # Handle matched samples with matched normal BAMs
            input_bams = snakemake.input[0:len(snakemake.input) - 3]
            input_bam = input_bams[:int(len(input_bams) / 2)]
            input_bai = input_bams[int(len(input_bams)/2):]
            
            inputs = snakemake.input[-3:len(snakemake.input)]
            batch_options = snakemake.params[4]

            # Print run info for logging
            print(f"Setting up batch scripts using the following inputs:\nBam files:\t{input_bam}\nBai files:\t{input_bai}\nParameters:\t{snakemake.params[6]}\nBatch options:\t{batch_options}")

            input_maf = open(inputs[0], "r")

            # Skip if no variants in outfile
            line_count = 0
            for line in input_maf:
                line_count += 1
                if line_count > 1:
                    break
            if line_count < 2:
                input_maf.close()
                touch_outputs(
                    output_dir = snakemake.params[0],
                    seq_type = snakemake.wildcards["seq_type"],
                    genome_build = snakemake.wildcards["genome_build"],
                    presets = snakemake.params[6],
                    tumour_id = snakemake.wildcards["tumour_id"],
                    suffix = snakemake.params[5],
                    finished_file = snakemake.output[0]
                )
                exit()

            # Return to top of MAF
            input_maf.seek(0)

            # Read MAF file and create dataframe
            regions = get_regions_df(
                input_maf,
                padding=batch_options["padding"]
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
                suffix = snakemake.params[5],
                igv_presets = snakemake.params[6],
                igv_options = batch_options["igv_options"],
                max_height = batch_options["max_height"],
                sleep_timer = batch_options["sleep_timer"]
            )

            touch_output = open(snakemake.output[0], "w")
            touch_output.close()
        
        except Exception as e:
            logging.error(e, exc_info=1)
            raise

def touch_outputs(output_dir, seq_type, genome_build, presets, tumour_id, suffix, finished_file):
    tumour_suffix = tumour_id + suffix + ".batch"
    for preset in presets:
        os.makedirs(os.path.join(output_dir, "--".join([seq_type, genome_build]), preset), exist_ok = True)
        merged_batch = os.path.join(output_dir, "merged_batch_scripts", "--".join([seq_type, genome_build]), preset, tumour_suffix)
        merged_file = open(merged_batch, "w")
        merged_file.close()
    touch_finished = open(finished_file, "w")
    touch_finished.close()

def get_regions_df(input_maf, padding):
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
        "padding": padding,
        "pair_status": maf.pair_status
        }
    )

    return regions_df

def output_lines(lines, batch_output):
    output = open(batch_output, "w")
    lines.append("")
    text = "\n".join(lines)
    output.write(text)
    output.close()

def generate_igv_batch_per_row(sleep_interval, preset, options, coordinates, directory, child_dir, seq_build, chrom_directory, snapshot_filename):
    lines = []

    lines.append(f"goto {coordinates}")

    snapshot_regions_dir = os.path.join(directory, seq_build, child_dir, preset, chrom_directory, "")

    lines.append(f"snapshotDirectory {snapshot_regions_dir}")
    for igv_option in options[preset]:
        lines.append(igv_option)
    lines.append(f"setSleepInterval {sleep_interval}")
    lines.append("collapse")
    lines.append(f"snapshot {snapshot_filename}")
    
    return lines

def generate_igv_batch_header(bam, index, max_height, genome_build):
    lines = []

    genome_build = genome_build.replace("grch37","hg19")

    assert len(bam) == len(index), "Error while generating batch script: number of .bam files and .bai files are not equal"

    for i in range(0,len(bam)):
        lines.append(f"load {bam[i]} index={index[i]}")

    lines.append(f"maxPanelHeight {max_height}")
    lines.append(f"genome {genome_build}")

    return lines

def generate_igv_batches(regions, bam, bai, output_dir, snapshot_dir, genome_build, seq_type, suffix, igv_presets, igv_options, max_height, sleep_timer=2000):
    for preset in igv_presets:
        for _, row in regions.iterrows():
            all_lines = []

            merged_batch_suffix = row.sample_id + suffix + ".batch"

            header = generate_igv_batch_header(bam=bam, index=bai, max_height=max_height, genome_build=genome_build)
            all_lines.extend(header)

            if row.pair_status == "matched":
                child_directory = "tumour_normal_pair"
            elif row.pair_status == "unmatched":
                child_directory = "tumour_only"

            seq_type_build = f"{seq_type}--{genome_build}"
            chrom_dir = row.chromosome

            filename = []
            filename.append(row.region),
            filename.append(row.region_name)
            filename.append(row.sample_id)

            batch_filename = "--".join(filename) + suffix + ".batch"
            filename = "--".join(filename) + suffix + ".png"

            lines = generate_igv_batch_per_row(
                sleep_interval = sleep_timer,
                preset = preset,
                options = igv_options,
                coordinates = row.snapshot_coordinates,
                directory = snapshot_dir,
                child_dir = child_directory,
                seq_build = seq_type_build,
                chrom_directory = chrom_dir,
                snapshot_filename = filename
            )

            all_lines.extend(lines)

            # Make subdirectories if necessary because snakemake won't make them since rule is a checkpoint
            os.makedirs(os.path.join(output_dir, "single_batch_scripts", seq_type_build, preset), exist_ok=True)

            batch_file_path = os.path.join(output_dir, "single_batch_scripts", seq_type_build, preset, batch_filename)

            output_lines(all_lines, batch_file_path)
        
        os.makedirs(os.path.join(output_dir, "merged_batch_scripts", seq_type_build, preset), exist_ok=True)

        merged_preset_path = os.path.join(output_dir, "merged_batch_scripts", seq_type_build, preset, merged_batch_suffix)

        merged_preset_touch = open(merged_preset_path, "w")
        merged_preset_touch.close()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[1],
        filemode='w'
    )
    main()
