#!/usr/bin/env python3

##### MODULES #####

import os
from os.path import join
from functools import reduce
from itertools import product
from collections import defaultdict
from glob import glob
import pandas as pd


##### UTILITIES #####

def ifelse(one, two):
    if one is None:
        return two
    return one


def symlink(source, target):
    target_dir = os.path.dirname(target)
    source_rel = os.path.relpath(source, target_dir)
    os.symlink(source_rel, target)


def get_from_dict(dictionary, access_list):
    """Access nested index in dictionary"""
    return reduce(dict.get, access_list, dictionary)


##### FILE SEARCHING #####

def locate_file(directory, file_ext, *patterns):
    """
    Locate a file with a given file extension in a given directory.
    The search is done recursively. Optionally filter for file paths 
    that contain the given patterns.
    """
    directory = directory.rstrip("/")
    file_ext  = file_ext.strip(".")
    patterns_fmt = "*" + "*, *".join(patterns) + "*" if patterns else None
    matches = []
    for root, subdirs, files in os.walk(directory):
        for subdir in subdirs:
            glob_pattern = f"{root}/{subdir}/*.{file_ext}"
            matches_here = glob(glob_pattern)
            matches.extend(matches_here)
    for pattern in patterns:
        matches = [x for x in matches if pattern in x]
    matches = [x for x in matches if "archive" not in x]
    if len(matches) == 0:
        msg = (f"No '.{file_ext}' file was found in '{directory}/' "
               f"with the following patterns: {patterns_fmt}")
        raise FileNotFoundError(msg)
    elif len(matches) > 1:
        msg = (f"More than one '.{file_ext}' file was found in '{directory}/' "
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
        A dictionary with the samples table under the "all" key.
    """
    # Load samples table into "all" key of dictionary
    samples = dict()
    samples["all"] = pd.read_table(file_path)
    # If a mapper function was provided, use it to rename columns
    if mapper:
        samples["all"].rename(columns=mapper, inplace=True)
    # If maps were provided, use them to rename columns
    if maps:
        maps_rev = {v: k for k, v in maps.items()}
        samples["all"].rename(columns=maps_rev, inplace=True)
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


def group_samples(samples, *subgroups):
    samples_dict = {}
    # Expect at least one subgroup
    if not subgroups:
        raise ValueError("Need to provide at least one subgroup.")
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
        sample_id = row["sample_id"]
        if sample_id in parent:
            raise ValueError(f"Sample ID ({sample_id}) not unique in this "
                             f"series of subgroups ({values}).")
        parent.append(sample_id)
    return samples_dict


def generate_pairs(samples):
    # Organize samples by patient and tissue status (tumour vs. normal)
    patients = group_samples(samples, "patient_id", "tissue_status")
    # Find every possible tumour-normal pair for each patient
    all_pairs = {
        "tumour": [],
        "normal": []
    }
    for patient, samples in patients.items():
        tumour_ids = samples.get("Tumour", [])
        normal_ids = samples.get("Normal", [])
        pairs = product(tumour_ids, normal_ids)
        for t, n in pairs:
            all_pairs["tumour"].append(t)
            all_pairs["normal"].append(n)
    return all_pairs


##### MODULE SETUP #####

def setup_module(config, name, version):
    
    # Get configuration for the given module
    mconfig = config['modules'][name]
    
    # Configure output directory if not specified
    if mconfig["output_dir"] is None:
        mconfig["output_dir"] = join(config["modules"]["output_dir"],
                                     f"{name}_module-{version}")
    
    # Figure out path to module sub-directory in repository
    msubdir = join(config['modules']["repository"], "modules", name, version)
    for env_name, env_vals in mconfig["conda_envs"].items():
        if env_vals["path"] is None:
            mconfig["conda_envs"][env_name] = join(msubdir, env_vals["yaml"])
        else:
            mconfig["conda_envs"][env_name] = env_vals["path"]
    
    # Return module-specific configuration
    return mconfig
