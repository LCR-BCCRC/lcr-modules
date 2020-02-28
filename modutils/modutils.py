#!+usr/bin/env python3

##### MODULES #####

import os
from os.path import isdir, join, realpath, relpath, split, normpath
from textwrap import dedent
from functools import reduce
from itertools import product
from collections import defaultdict, namedtuple
from glob import glob
import pandas as pd
from yaml import dump
from snakemake.logging import logger


##### UTILITIES #####

def symlink(src, dest):
    """
    Create relative symlink to file from any working directory.
    """
    if isdir(dest):
        dest_dir  = dest
        dest_file = split(src)[1]
    else:
        dest_dir, dest_file = split(dest)
    src_rel = relpath(src, dest_dir)
    dest = join(dest_dir, dest_file)
    os.symlink(src_rel, dest)


def get_from_dict(dictionary, access_list):
    """Access nested index in dictionary"""
    return reduce(dict.get, access_list, dictionary)


def collapse(text):
    """Collapse a triple-quoted multi-line string to one line"""
    lines = text.strip().split("\n")
    lines_dedented = [x.strip(" ") for x in lines]
    return " ".join(lines_dedented)


##### FILE SEARCHING #####

def locate_file(directory, file_ext, *patterns):
    """
    Locate a file with a given file extension in a given directory.
    The search is done recursively. Optionally filter for file paths 
    that contain the given patterns.
    """
    file_ext  = file_ext.strip(".")
    patterns_fmt = "*" + "*, *".join(patterns) + "*" if patterns else None
    matches = []
    for root, subdirs, files in os.walk(directory):
        for subdir in subdirs:
            glob_pattern = join(root, subdir, f"*.{file_ext}")
            matches_here = glob(glob_pattern)
            matches.extend(matches_here)
    for pattern in patterns:
        matches = [x for x in matches if pattern in x]
    matches = [x for x in matches if "archive" not in x]
    if len(matches) == 0:
        msg = (f"No '.{file_ext}' file was found in '{directory}' "
               f"with the following patterns: {patterns_fmt}")
        raise FileNotFoundError(msg)
    elif len(matches) > 1:
        msg = (f"More than one '.{file_ext}' file was found in '{directory}' "
               f"with the following patterns: {patterns_fmt}")
        raise FileNotFoundError(msg)
    else:
        return matches[0]


def locate_bam(*patterns):
    bam = locate_file("data", "bam", *patterns)
    return bam


def locate_bams_generator(seq_type):
    def locate_bams_custom(wildcards):
        keys = wildcards.keys()
        bams = dict()
        if "sample_id" in keys:
            sample_bam = locate_bam(wildcards.sample_id, seq_type)
            bams["sample_bam"] = sample_bam
        else:
            tumour_bam = locate_bam(wildcards.tumour_id, seq_type)
            normal_bam = locate_bam(wildcards.normal_id, seq_type)
            bams["tumour_bam"] = tumour_bam
            bams["normal_bam"] = normal_bam
        return bams
    return locate_bams_custom


locate_genome_bams  = locate_bams_generator("genome")
locate_mrna_bams    = locate_bams_generator("mrna")
locate_mirna_bams   = locate_bams_generator("mirna")
locate_capture_bams = locate_bams_generator("capture")


##### SAMPLE PROCESSING #####

def load_samples(file_path, mapper=None, **maps):
    """
    Load samples TSV file and rename columns using a `mapper` function,
    a series of `maps` pairs (after="before"), or both.
    
    Returns:
        A pandas data frame.
    """
    # Load samples table into "all" key of dictionary
    samples = pd.read_table(file_path)
    # If a mapper function was provided, use it to rename columns
    if mapper:
        samples.rename(columns=mapper, inplace=True)
    # If maps were provided, use them to rename columns
    if maps:
        maps_rev = {v: k for k, v in maps.items()}
        samples.rename(columns=maps_rev, inplace=True)
    return samples


def filter_samples(samples, **filters):
    """
    Filter a samples table based on the provided `filters`.
    The key of a filter refers to the column to filter on. 
    The value refers to the filter value. THe value can
    either be a single value or a list of values. 
    
    Returns:
        A filtered table. 
    """
    for column, value in filters.items():
        if isinstance(value, list):
            samples = samples[samples[column].isin(value)]
        else:
            samples = samples[samples[column] == value]
    return samples


def group_samples(samples, 
                  subgroups=["seq_type", "patient_id"], 
                  features=["sample_id"]):
    # Setup
    samples_dict = dict()
    Sample = namedtuple("Sample", features)
    # Expect at least one subgroup
    if len(subgroups) == 0:
        raise ValueError("Need to provide at least one subgroup.")
    # Expect at least one feature
    if len(features) == 0:
        raise ValueError("Need to provide at least one feature.")
    # Iterate over each row
    for index, row in samples.iterrows():
        values = []
        # Initialize intermediate subgroups with dictionaries
        parent = samples_dict
        for subgroup in subgroups[:-1]:
            value = row[subgroup]
            if value not in parent:
                parent[value] = {}
            values.append(value)
            parent = get_from_dict(samples_dict, values)
        # Initialize "terminal" subgroups with sets
        value = row[subgroups[-1]]
        values.append(value)
        if value not in parent:
            parent[value] = list()
        # Add sample ID to the "terminal" subgroup
        parent = parent[value]
        sample = Sample(*row[features])
        if sample in parent:
            raise ValueError(f"`{sample}` not unique for this "
                             f"nested set of subgroups ({values}).")
        parent.append(sample)
    return samples_dict


def generate_runs_for_patient(samples, return_paired=True, 
                              return_unpaired=False, 
                              return_paired_as_unpaired=False):
    runs = defaultdict(list)
    tumour_samples = samples.get("Tumour", [])
    normal_samples = samples.get("Normal", [None])
    # If return_paired_as_unpaired is True, then return_unpaired should be True
    return_unpaired = return_unpaired or return_paired_as_unpaired
    # Add an unpaired normal is there isn't one
    if return_paired_as_unpaired and None not in normal_samples:
        normal_samples.append(None)
    pairs = product(tumour_samples, normal_samples)
    for tumour, normal in pairs:
        # Don't return paired runs if normal is present
        if not return_paired and normal is not None:
            continue
        # Don't return unpaired runs if normal is absent
        if not return_unpaired and normal is None:
            continue
        # Compile features
        tumour = tumour._asdict()
        if normal is None:
            normal = defaultdict(lambda: None)
        else:
            normal = normal._asdict()
        for field in tumour.keys():
            runs["tumour_" + field].append(tumour[field])
            runs["normal_" + field].append(normal[field])
    return dict(runs)


def combine_lists(dictionary):
    combined = defaultdict(list)
    for d in dictionary.values():
        for k, v in d.items():
            combined[k].extend(v)
    combined_plain = dict(combined)
    combined_df = pd.DataFrame(combined_plain)
    return combined_df


def walk_through_dict(dictionary, end_fn=print, end_depth=None,
                      trace=None, result=None, **kwargs):
    # Define default values
    if end_depth is None:
        end_depth = float('inf')
    if trace is None:
        trace = tuple()
    if result is None:
        result = dict()
    # Begin recursion
    if end_depth <= 0:
        return end_fn(dictionary, **kwargs)
    for k, v in dictionary.items():
        current_trace = trace + (k,)
        if isinstance(v, dict) and len(current_trace) < end_depth:
            result[k] = {}
            walk_through_dict(v, end_fn, end_depth, current_trace, 
                              result[k], **kwargs)
        else:
            result[k] = end_fn(v, **kwargs)
        if len(result[k]) == 0:
            del result[k]
    return result


def generate_runs(samples,
                  subgroups=["seq_type"],
                  features=["sample_id"],
                  **kwargs):
    # Copy subgroups to avoid using the same mutable object in every call
    subgroups = subgroups.copy()
    # Ensure that patient_id and tissue_status are the last two subgroups
    subgroups.extend(["patient_id", "tissue_status"])
    # Organize samples by patient and tissue status (tumour vs. normal)
    patients = group_samples(samples, subgroups, features)
    # Find every possible tumour-normal pair for each patient
    end_depth = len(subgroups) - 1
    runs_by_patient = walk_through_dict(patients, generate_runs_for_patient, 
                                        end_depth, **kwargs)
    runs = walk_through_dict(runs_by_patient, combine_lists, end_depth - 1)
    return runs


##### MODULE SETUP/CLEANUP #####

def setup_module(config, name, version):
    
    # Get configuration for the given module
    mconfig = config["modules"][name]
    
    # Find repository and module directories
    repodir = normpath(config["modules"]["_shared"]["repository"])
    msubdir = join(repodir, "modules", name, version)
    
    # Ensure that common module sub-fields are present
    subfields = ["inputs", "outputs", "dirs", "conda_envs", 
                 "options", "threads", "memory"]
    for subfield in subfields:
        if subfield not in mconfig:
            logger.warning(f"Subfield '{subfield}' is missing from config "
                           f"for module {name}")
            mconfig[subfield] = {}
    
    # Update file paths with "{REPODIR}" to point to the repository directory
    for pairs in mconfig.values():
        for k, v in pairs.items():
            if isinstance(v, str) and "{REPODIR}" in v:
                pairs[k] = v.replace("{REPODIR}", repodir)
    
    # Configure output directory if not specified and create it
    if mconfig["dirs"].get("_parent") is None:
        root_output_dir = config["modules"]["_shared"].get("root_output_dir")
        root_output_dir = root_output_dir or "modules"
        output_dir = join(root_output_dir,  f"{name}-{version}")
        mconfig["dirs"]["_parent"] = output_dir
    os.makedirs(mconfig["dirs"]["_parent"], exist_ok=True)
    
    # Update paths to conda environments to be relative to the module directory
    for env_name, env_val in mconfig["conda_envs"].items():
        if env_val is not None:
            mconfig["conda_envs"][env_name] = relpath(env_val, msubdir)
    
    # Normalize all paths
    for pairs in mconfig.values():
        for k, v in pairs.items():
            if isinstance(v, str) and "{REPODIR}" in v:
                pairs[k] = v.replace("{REPODIR}", repodir)
    
    # Return module-specific configuration
    return mconfig


def setup_subdirs(module_config, *subdirs):
    if "_parent" in subdirs:
        raise ValueError("You cannot have a sub-directory called '_parent'. "
                         "Consider using a more specific term.")
    subdirs = ["input", *subdirs, "output"]
    numbers = [ f"{x:02}" for x in (*range(0, len(subdirs) - 1), 99) ]
    for num, subdir in zip(numbers, subdirs):
        subdir_full = join(module_config["dirs"]["_parent"], f"{num}-{subdir}")
        module_config["dirs"][subdir] = subdir_full
        os.makedirs(subdir_full, exist_ok=True)
    return module_config


def cleanup_module(module_config):
    # Define useful variables
    parent_dir = module_config["dirs"]["_parent"]
    # Define fields to be output as TSV files
    tsv_fields = {"samples": None, "runs": None}
    for field in tsv_fields.keys():
        tsv_fields[field] = module_config.pop(field)
        output_file = join(parent_dir, f"{field}.tsv")
        tsv_fields[field].to_csv(output_file, sep="\t", index=False)
    # Output current configuration for future reference
    config_file = join(parent_dir, "config.yaml")
    with open(config_file, "w") as config_file_handler:
        dump(module_config, config_file_handler)
    # Add back the TSV fields
    for field in tsv_fields.keys():
        module_config[field] = tsv_fields[field]
    # Return nothing
    return None
