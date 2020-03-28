#!/usr/bin/env python3

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
from snakemake.utils import min_version, validate


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


##### SNAKEMAKE INPUT/PARAM FUNCTIONS #####

def make_seqtype_specific(param_config):
    """
    Return parameter function to be used by Snakemake to
    retrieve the correct parameters based on the seq_type.
    """
    
    def run_seqtype_param(wildcards):
        param = [
            param_config.get("_base", ""),
            param_config.get(wildcards.seq_type, "")
        ]
        # Stich non-null parameters
        param_str = " ".join(x for x in param if x != "")
        return param_str
    
    return run_seqtype_param


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

LOWERCASE_COLS = ("tissue_status", "seq_type", "ff_or_ffpe")

def load_samples(file_path, lowercase_cols=LOWERCASE_COLS, mapper=None, **maps):
    """
    Load samples TSV file and rename columns using a `mapper` function,
    a series of `maps` pairs (after="before"), or both.
    
    As a convenience feature, this function will also force a few
    columns to lower case to avoid issues downstream. The default
    set of columns to be forced to lower case are defined in:
    
        LOWERCASE_COLS
    
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
    # Force a few important columns to be lowercase
    for col in lowercase_cols:
        if col not in samples.columns:
            logger.warn(f"The `{col}` column does not exist "
                        "in the samples data frame.")
            continue
        samples[col] = samples[col].str.lower()
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
        if isinstance(value, str):
            value = [value]
        samples = samples[samples[column].isin(value)]
    return samples


def group_samples(samples, sample_class, subgroups=["seq_type"]):
    # Setup
    samples_dict = dict()
    # Expect at least one subgroup
    if len(subgroups) == 0:
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
        sample = sample_class(*row)
        if sample in parent:
            raise ValueError(f"`{sample}` not unique for this "
                             f"nested set of subgroups ({values}).")
        parent.append(sample)
    return samples_dict


def generate_runs_for_patient(samples, run_unpaired_tumours_with=None,
                              unmatched_normal=None, run_paired_tumour=True,
                              run_paired_tumours_as_unpaired=False):
    """Generate a run for every tumour with and/or without a paired normal.
    
    Note that 'unpaired tumours' in the argument names and documentation
    refers to tumours without a matched normal sample. 
    
    Args:
        samples: A dict containing lists of tumour and normal sample IDs (str)
            for a given patient (in any order).
        run_unpaired_tumours_with: A str specifying what to pair with unpaired
            tumours. Value must be one of the following values: 
                None, 'no_normal', or 'unmatched_normal'
            This must be set if run_paired_tumours_as_unpaired is True.
        unmatched_normal: A NamedTuple representing the normal sample to be
            used with unpaired tumours when `run_unpaired_tumours_with` is set
            to 'unmatched_normal'.
        run_paired_tumour: A bool specifying whether to run paired tumours.
            This is useful for naturally unpaired analyses (e.g., for RNA-seq).
        run_paired_tumours_as_unpaired: A bool specifying whether paired 
            tumours should also be run as unpaired (i.e., separate from their
            matched normal sample). This is useful for benchmarking purposes.
    
    Returns:
        A dict containing lists of sample features (extracted from the samples
        data frame) prefixed with 'tumour_' and 'normal_' for all tumours for
        the given patient. Depending on the argument values, tumour-normal
        pairs may not be matching, and normal samples may not be included.
        A column called 'pair_status' specifies whether a tumour is paired with
        a matched normal sample.
    """

    run_unpaired_tumours_with_options = (None, "no_normal", "unmatched_normal")
    if run_unpaired_tumours_with not in run_unpaired_tumours_with_options:
        msg = ("`run_unpaired_tumours_with` must be one of the values below "
               f"(not `{run_unpaired_tumours_with}`): \n"
               f"{run_unpaired_tumours_with_options}")
        raise ValueError(msg)
    
    run_unpaired_tumour = run_unpaired_tumours_with is not None
    
    if run_paired_tumours_as_unpaired and not run_unpaired_tumour:
        msg = ("`run_paired_tumours_as_unpaired` was True whereas "
               "`run_unpaired_tumours_with` was None. Please set "
               "`run_unpaired_tumours_with` to 'unmatched_normal' "
               "or 'no_normal'.")
        raise ValueError(msg)
    
    if (run_unpaired_tumours_with == "unmatched_normal" 
            and unmatched_normal is None):
        msg = ("`run_unpaired_tumours_with` was set to 'unmatched_normal' "
               "whereas `unmatched_normal` was None. For each seq_type, "
               "provide an unmatched normal sample ID in the _shared section "
               "of the modules configuration under `unmatched_normal_id`.")
        raise ValueError(msg)
    
    # Retrieve tumour and normal samples
    runs = defaultdict(list)
    tumour_samples = samples.get("tumour", []) + samples.get("tumor", [])
    normal_samples = samples.get("normal", [None])
    
    # Add an unpaired normal is there isn't one
    if run_paired_tumours_as_unpaired and None not in normal_samples:
        normal_samples.append(None)
    
    for tumour, normal in product(tumour_samples, normal_samples):
        # Check for paired samples
        paired = normal is not None
        if paired and run_paired_tumour is False:
            continue
        # Check for unpaired samples
        unpaired = normal is None
        if unpaired and run_unpaired_tumour is False:
            continue
        # Compile features
        tumour = tumour._asdict()
        if normal is None and run_unpaired_tumours_with == "unmatched_normal":
            normal = unmatched_normal._asdict()
            runs["pair_status"].append("unmatched")
        elif normal is None and run_unpaired_tumours_with == "no_normal":
            normal = { key: None for key in tumour.keys() }
            runs["pair_status"].append("no_normal")
        else:
            normal = normal._asdict()
            runs["pair_status"].append("matched")
        for field in tumour.keys():
            runs["tumour_" + field].append(tumour[field])
            runs["normal_" + field].append(normal[field])
    
    return dict(runs)


def generate_runs_for_patient_and_seqtype(samples, pairing_config, **kwargs):
    """Run generate_runs_for_patient based on the current seq_type."""
    
    seq_type_set = set()
    for samples_list in samples.values():
        seq_type_set.update(s.seq_type for s in samples_list)
    
    if len(seq_type_set) != 1:
        msg = ("This function is only meant to be run on groups of samples "
               "for a given patient and a given sequencing data type. The "
               "current group of samples has the following seq_types: \n"
               f"{seq_type}")
        raise ValueError(msg)
    
    seq_type = seq_type_set.pop()
    pconfig = pairing_config[seq_type]
    unmatched_normals = kwargs.pop("unmatched_normals")
    
    if pconfig["run_unpaired_tumours_with"] == "unmatched_normal":
        unmatched_normal  = unmatched_normals[seq_type]
        runs = generate_runs_for_patient(samples, **pconfig, **kwargs,
                                        unmatched_normal=unmatched_normal)
    else:
        runs = generate_runs_for_patient(samples, **pconfig, **kwargs)
    
    return runs


def combine_lists(dictionary):
    combined = defaultdict(list)
    for d in dictionary.values():
        for k, v in d.items():
            combined[k].extend(v)
    combined = dict(combined)
    combined = pd.DataFrame(combined)
    return combined


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


def generate_runs(samples, sample_class, subgroups=[], **kwargs):
    
    # Copy subgroups to avoid using the same mutable object in every call
    subgroups = subgroups.copy()
    
    # Ensure that patient_id and tissue_status are the last two subgroups
    subgroups.extend(["patient_id", "tissue_status"])
    
    # Organize samples by patient and tissue status (tumour vs. normal)
    patients = group_samples(samples, sample_class, subgroups)
    
    # Find every possible tumour-normal pair for each patient
    end_depth = len(subgroups) - 1
    runs = walk_through_dict(patients, generate_runs_for_patient_and_seqtype, 
                             end_depth, **kwargs)
    while end_depth > 0:
        runs = walk_through_dict(runs, combine_lists, end_depth - 1)
        end_depth -= 1
    
    # Warn if runs have duplicates
    if any(runs.duplicated()):
        logger.warn("Duplicate runs exist. This probably shouldn't happen.")
    
    return runs


##### MODULE SETUP/CLEANUP #####

def setup_module(config, name, version, subdirs, **kwargs):
    
    # Ensure minimum version of Snakemake
    min_version("5.0.0")
    
    # Get configuration for the given module and create samples shorthand
    mconfig = config["modules"][name]
    sconfig = config["modules"]["_shared"]
    msamples = mconfig["samples"]
    
    # Find repository and module directories
    repodir = normpath(sconfig["repository"])
    msubdir = join(repodir, "modules", name, version)
    
    # Ensure that common module sub-fields are present
    subfields = ["inputs", "dirs", "conda_envs", 
                 "options", "threads", "memory"]
    for subfield in subfields:
        if subfield not in mconfig:
            logger.warning(f"Subfield '{subfield}' is missing from config "
                           f"for module {name}")
            mconfig[subfield] = dict()
    
    # Update file paths with "{REPODIR}" to point to the repository directory
    for pairs in mconfig.values():
        for k, v in pairs.items():
            if isinstance(v, str) and "{REPODIR}" in v:
                pairs[k] = v.replace("{REPODIR}", repodir)
    
    # Validation samples data frame
    schemas_dir = join(msubdir, "schemas")
    schemas = os.listdir(schemas_dir)
    for schema in schemas:
        validate(msamples, schema=join(schemas_dir, schema))
    
    # Configure output directory if not specified and create it
    if mconfig["dirs"].get("_parent") is None:
        root_output_dir = sconfig.get("root_output_dir")
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
    
    # Setup sub-directories
    mconfig = setup_subdirs(mconfig, subdirs)
    
    # Setup for unpaired tumour samples
    Sample = namedtuple("Sample", msamples.columns.tolist())
    unmatched_normals = None
    unmatched_normal_ids = sconfig.get("unmatched_normal_ids")
    if unmatched_normal_ids is not None:
        unmatched_normals = dict()
        for seq_type, normal_id in unmatched_normal_ids.items():
            normal_row = msamples[msamples.sample_id == normal_id]
            num_matches = len(normal_row)
            if num_matches != 1:
                msg = (f"There are {num_matches} samples matching the normal "
                       f"ID ({normal_id}) rather than expected (1)")
                raise ValueError(msg)
            unmatched_normals[seq_type] = Sample(*normal_row.squeeze())
    
    # Generate runs
    runs = generate_runs(msamples, sample_class=Sample, subgroups=["seq_type"],
                         pairing_config=mconfig["pairing_config"],
                         unmatched_normals=unmatched_normals, **kwargs)
    
    # Split runs based on pair_status
    mconfig["runs"] = runs
    mconfig["paired_runs"]   = runs[runs.pair_status != "no_normal"]
    mconfig["unpaired_runs"] = runs[runs.pair_status == "no_normal"]
    
    # Return module-specific configuration
    return mconfig


def setup_subdirs(module_config, subdirs):
    if "_parent" in subdirs:
        raise ValueError("You cannot have a sub-directory called '_parent'. "
                         "Consider using a more specific term.")
    subdirs = ["inputs", *subdirs, "outputs"]
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
    tsv_fields = {"samples": None, "paired_runs": None, "unpaired_runs": None}
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
