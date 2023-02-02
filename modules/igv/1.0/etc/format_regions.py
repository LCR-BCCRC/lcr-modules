#!/usr/bin/env python

import os
import pandas as pd
import oncopipe as op

def format_clustl(clustl_regions):
    # Convert OncodriveCLUSTL cluster coordinates to BED format
    p_filter = CLUSTL_PARAMS["p_value"]
    score_filter = CLUSTL_PARAMS["score"]
    n_samples_filter = CLUSTL_PARAMS["n_samples"]

    for key, filter_value in {"P": p_filter, "SCORE": score_filter, "N_SAMPLES": n_samples_filter}.items():
        if filter_value is not None:
            if key != "P":
                clustl_regions = clustl_regions[clustl_regions[key] >= float(filter_value)]
            if key == "P":
                clustl_regions = clustl_regions[clustl_regions[key] <= float(filter_value)]

    # Reformat CLUSTL coordinates to handle clusters that cross introns (when CLUSTL concatenated mode is used)
    clustl_regions = clustl_regions.assign(COORDINATES = clustl_regions.COORDINATES.str.split(";")).explode("COORDINATES")
    clustl_regions["COORDINATES"] = clustl_regions.apply(
        lambda x: list(
            range(
                int(str(x["COORDINATES"]).split(",")[0]), int(str(x["COORDINATES"]).split(",")[1]) + 1
            )
        )
        if str(x["COORDINATES"]).split(",")[0] != str(x["COORDINATES"]).split(",")[1] else int(str(x["COORDINATES"]).split(",")[0]),
        axis = 1
    )
    clustl_regions = clustl_regions.explode("COORDINATES")

    # Create columnsn required for BED format
    chr_str = "chr" + clustl_regions["CHROMOSOME"].map(str)
    clustl_reformatted = pd.DataFrame(
        {
            "chrom": chr_str,
            "start": clustl_regions["COORDINATES"],
            "end": clustl_regions["COORDINATES"]
        }
    )
    return clustl_reformatted

def format_maf(regions):
    # If the regions format is a MAF, don't need to reformat for liftover
    return regions

def format_regions(regions, regions_format):
    format_functions = {
        "maf": format_maf,
        "oncodriveclustl": format_clustl,
        "genomic_pos": format_genomic_pos
    }

    return format_functions[regions_format](regions)

regions_file = snakemake.input[0]
regions_format = snakemake.params[0]

output_file = snakemake.output[0]

if regions_format == "oncodriveclustl":
    CLUSTL_PARAMS = snakemake.params[1]

# Read regions into dataframe
regions_df = pd.read_table(regions_file, comment="#", sep="\t")

regions_formatted = format_regions(regions_df, regions_format)

regions_formatted.to_csv(output_file, sep="\t", index=False)