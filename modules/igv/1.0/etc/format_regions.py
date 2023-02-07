#!/usr/bin/env python

import os
import pandas as pd
import oncopipe as op

def format_mutation_id(mutation_id):
    ## Modify dataframe to handle NA values in columns
    #mutation_id = mutation_id.fillna(0)
#
    ## Filter dataframe based on config option
    #max_distinct_genome_cohorts = MUTATION_ID_PARAMS["max_distinct_genome_cohorts"]
    #min_total_genome_samples = MUTATION_ID_PARAMS["min_total_genome_samples"]
    #
    #max_distinct_capture_cohorts = MUTATION_ID_PARAMS["max_distinct_capture_cohorts"]
    #min_total_capture_samples = MUTATION_ID_PARAMS["min_total_capture_samples"]
#
    #for key, filter_value in {
    #    "distinct_genome_cohorts": max_distinct_genome_cohorts, 
    #    "total_genome": min_total_genome_samples, 
    #    "distinct_capture_cohorts": max_distinct_capture_cohorts,
    #    "total_capture": min_total_capture_samples
    #}.items():
    #    if filter_value is not None:
    #        if key in ["distinct_genome_cohorts", "distinct_capture_cohorts"]:
    #            mutation_id = mutation_id[mutation_id[key] <= float(filter_value)]
    #        else:
    #            mutation_id = mutation_id[mutation_id[key] >= float(filter_value)]

    # Create columns required for liftover in BED format
    genomic_pos_col = f"mutation_id_{REGIONS_BUILD}"

    for col, idx in {"chr_std": 0, "start": 1, "end": 2}.items():
        mutation_id[col] = mutation_id.apply(lambda x: str(x[genomic_pos_col]).split(":")[idx].replace("chr",""), axis=1)

    mutation_id_reformatted = pd.DataFrame(
        {
            "chrom": "chr" + mutation_id["chr_std"],
            "start": mutation_id["start"],
            "end": mutation_id["end"]
        }
    )

    # Remove duplicate rows
    mutation_id_reformatted = mutation_id_reformatted.drop_duplicates(keep='first')

    return mutation_id_reformatted

def format_hotmaps(hotmaps_regions):
    # Convert HotMAPS coordinates to BED format

    hotmaps_regions["chr_std"] = hotmaps_regions.apply(lambda x: str(x["Chromosome"]).replace("chr",""), axis=1)
    chr_std = "chr" + hotmaps_regions["Chromosome"].map(str)

    hotmaps_reformatted = pd.DataFrame(
        {
            "chrom": chr_std,
            "start": hotmaps_regions["Start_Position"],
            "end": hotmaps_regions["Start_Position"]
        }
    )
    return hotmaps_reformatted

def format_clustl(clustl_regions):
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
    
    # Convert OncodriveCLUSTL cluster coordinates to BED format
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
        "hotmaps": format_hotmaps,
        "mutation_id": format_mutation_id
    }

    return format_functions[regions_format](regions)

regions_file = snakemake.input[0]
regions_format = snakemake.params[0]

output_file = snakemake.output[0]

if regions_format == "oncodriveclustl":
    CLUSTL_PARAMS = snakemake.params[1]

if regions_format == "mutation_id":
    REGIONS_BUILD = snakemake.params[2]
    REGIONS_BUILD = REGIONS_BUILD.lower()

# Read regions into dataframe
regions_df = pd.read_table(regions_file, comment="#", sep="\t")

# Reformat for liftover based on regions format
regions_formatted = format_regions(regions_df, regions_format)

# Output regions file
regions_formatted.to_csv(output_file, sep="\t", index=False)