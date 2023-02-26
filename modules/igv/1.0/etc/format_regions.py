#!/usr/bin/env python

import os
import pandas as pd
import oncopipe as op
import vcf
import shutil

def format_mutation_id(mutation_id):
    # Read regions into dataframe
    mutation_id = pd.read_table(mutation_id, comment="#", sep="\t")

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
    # Read regions into dataframe
    hotmaps_regions = pd.read_table(hotmaps_regions, comment="#", sep="\t")

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
    # Read regions into dataframe
    clustl_regions = pd.read_table(clustl_regions, comment="#", sep="\t")

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

def format_vcf(regions):
    # Load VCF file
    vcf_reader = vcf.Reader(open(regions, "rb"))

    # Convert VCF records to BED format
    chroms = []
    pos = []
    events_seen = set()

    for record in vcf_reader:
        if len(record.FILTER) > 0:
            continue
        
        # Skip SVs with ID matching previous record
        if record.ID in events_seen:
            continue

        chromosome = "chr" + str(record.CHROM).replace("chr","")
        position = record.POS

        chroms.append(chromosome)
        pos.append(position)

        if record.is_sv and "END" in record.INFO:
            # Add end position of SV to regions of interest
            end = record.INFO["END"][0]

            chroms.append(chromosome)
            pos.append(end)

        if record.is_sv and record.INFO["SVTYPE"] == "BND":
            # Add end position of SV to regions of interest
            chromosome = "chr" + str(record.ALT[0].chr).replace("chr","")
            position = record.ALT[0].pos

            chroms.append(chromosome)
            pos.append(position)

            # To skip mate event in VCF file
            events_seen.add(record.INFO["MATEID"])
    
    vcf_reformatted = pd.DataFrame(
        {
            "chrom": chroms,
            "start": pos,
            "end": pos
        }
    )

    return vcf_reformatted

def format_regions(regions, regions_format):
    format_functions = {
        "oncodriveclustl": format_clustl,
        "hotmaps": format_hotmaps,
        "mutation_id": format_mutation_id,
        "vcf": format_vcf,
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

if regions_format == "bed" or regions_format == "maf":
    # Do not need to reformat for liftover
    shutil.copy(regions_file, output_file)
    exit()

# Reformat for liftover based on regions format
regions_formatted = format_regions(regions_file, regions_format)

# Output regions file
regions_formatted.to_csv(output_file, sep="\t", index=False)
