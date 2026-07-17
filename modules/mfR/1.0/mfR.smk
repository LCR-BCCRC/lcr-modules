#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  <your name>
# Module Author:    <your name>
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")
# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe
current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")
# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mfR"]`
CFG = op.setup_module(
    name = "mfR",
    version = "1.0",
    subdirectories = ["inputs", "extract", "foci", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mfR_input_maf,
    _mfR_input_subsets,
    _mfR_aggregate,
    _mfR_output_tsv,
    _mfR_all,


##### SAMPLE_SET -> SAMPLES MAPPING #####


# Sample_set membership lives in an external TSV, not in CFG["runs"], so it has
# to be resolved once here at parse time (not deferred to a script) so that the
# per-sample_set x chromosome extraction rule below can list exactly which
# per-sample MAFs it needs as its `input:`.
#
# Pre-capture these — CFG is deleted by op.cleanup_module and would be out of
# scope when Snakemake evaluates lambdas at DAG build time.
import pandas as pd
from collections import defaultdict

_sid_col = CFG["options"]["sample_id_column"]
_ss_col  = CFG["options"]["sample_set_column"]

_sample_sets_df = pd.read_csv(CFG["inputs"]["sample_sets"], sep = "\t")
if not {_sid_col, _ss_col}.issubset(_sample_sets_df.columns):
    raise ValueError(
        f"sample_sets TSV ({CFG['inputs']['sample_sets']}) must contain columns "
        f"'{_sid_col}' and '{_ss_col}' (see options.sample_id_column / options.sample_set_column)"
    )

# sample_id -> sample_set (a sample belongs to exactly one sample_set)
_sample_to_set = dict(zip(_sample_sets_df[_sid_col], _sample_sets_df[_ss_col]))

# sample_set -> list of {tumour_id, normal_id, seq_type, genome_build, pair_status}
# for every sample in CFG["runs"] that belongs to that sample_set. Assumes sample
# ids in the sample_sets TSV match CFG["runs"]["tumour_sample_id"] (i.e. the
# project's main sample table uses the same sample id convention).
_mfR_sampleset_runs = defaultdict(list)
for _t, _n, _s, _g, _p in zip(
    CFG["runs"]["tumour_sample_id"],
    CFG["runs"]["normal_sample_id"],
    CFG["runs"]["tumour_seq_type"],
    CFG["runs"]["tumour_genome_build"],
    CFG["runs"]["pair_status"],
):
    _sample_set = _sample_to_set.get(_t)
    if _sample_set is not None:
        _mfR_sampleset_runs[_sample_set].append(
            {"tumour_id": _t, "normal_id": _n, "seq_type": _s, "genome_build": _g, "pair_status": _p}
        )
# (not deleting the loop variables here: CFG["runs"] can legitimately be empty
# for some project configs, which would leave them unassigned and make a
# `del` raise NameError instead of the clearer "No samples found for
# sample_set" error raised later in _mfR_runs_for_sampleset)


##### FUNCTIONS #####


def _mfR_gather_chroms(wildcards):
    """Gather the per-chromosome foci tables for one sample_set."""
    return expand(
        CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.foci.tsv",
        sample_set = wildcards.sample_set,
        chrom = CFG["chromosomes"]
    )


def _mfR_runs_for_sampleset(wildcards):
    runs = _mfR_sampleset_runs.get(wildcards.sample_set, [])
    if not runs:
        raise ValueError(
            f"No samples found for sample_set '{wildcards.sample_set}' in "
            f"{CFG['inputs']['sample_sets']}. Check that sample_set values in the "
            f"config match values in the sample_set column of that TSV."
        )
    return runs


def _mfR_input_mafs_for_sampleset(wildcards):
    runs = _mfR_runs_for_sampleset(wildcards)
    return expand(
        str(rules._mfR_input_maf.output.maf_gz),
        zip,
        seq_type     = [r["seq_type"] for r in runs],
        genome_build = [r["genome_build"] for r in runs],
        tumour_id    = [r["tumour_id"] for r in runs],
        normal_id    = [r["normal_id"] for r in runs],
        pair_status  = [r["pair_status"] for r in runs],
    )


def _mfR_input_tbis_for_sampleset(wildcards):
    return [m + ".tbi" for m in _mfR_input_mafs_for_sampleset(wildcards)]


##### RULES #####


# Symlinks each sample's bgzipped + tabix-indexed genome MAF (and its .tbi)
# into the module results directory (under '00-inputs/'). Symlinked once per
# unique sample regardless of how many sample_sets it belongs to.
rule _mfR_input_maf:
    input:
        maf_gz = CFG["inputs"]["sample_maf"]
    output:
        maf_gz = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf.gz",
        tbi    = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf.gz.tbi"
    run:
        op.absolute_symlink(input.maf_gz, output.maf_gz)
        op.absolute_symlink(input.maf_gz + ".tbi", output.tbi)


# Symlinks the sample-set membership table into the module results directory
rule _mfR_input_subsets:
    input:
        sample_sets = CFG["inputs"]["sample_sets"]
    output:
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)


# For one sample_set x chromosome: tabix-extract just that chromosome's rows
# from just that sample_set's per-sample MAFs (never touches other samples or
# other chromosomes), and drop coding Variant_Classification rows. Replaces
# the old prepare_maf.R (whole-master-MAF read) entirely, and replaces the
# whole-per-sample_set-MAF read that used to happen once per chromosome.
rule _mfR_extract_chrom:
    input:
        mafs = _mfR_input_mafs_for_sampleset,
        tbis = _mfR_input_tbis_for_sampleset
    output:
        maf = temp(CFG["dirs"]["extract"] + "{sample_set}/chromosomes/{chrom}.noncoding.maf")
    log:
        CFG["logs"]["extract"] + "{sample_set}/{chrom}.extract.log"
    conda:
        CFG["conda_envs"]["mfr"]
    container:
        CFG["container_envs"]["mfr"]
    threads:
        CFG["threads"]["extract_chrom"]
    resources:
        **CFG["resources"]["extract_chrom"]
    params:
        coding_classes = CFG["options"]["coding_variant_classifications"]
    script:
        "src/python/extract_chrom_maf.py"


# Cluster non-coding mutation positions into "foci", scattered by chromosome.
# This is the heavy step (one job per sample_set x chromosome). The input is
# already scoped to this sample_set x chromosome, so no further filtering of
# samples or chromosomes happens here — just the clustering math.
rule _mfR_cluster:
    input:
        maf = str(rules._mfR_extract_chrom.output.maf)
    output:
        tsv  = CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.foci.tsv",
        plot = CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.silhouette.pdf"
    log:
        log = CFG["logs"]["foci"] + "{sample_set}/{chrom}.cluster.log"
    conda:
        CFG["conda_envs"]["mfr"]
    container:
        CFG["container_envs"]["mfr"]
    threads:
        CFG["threads"]["cluster"]
    resources:
        # dist() is O(n^2) in unique positions; bump mem for dense chromosomes.
        **CFG["resources"]["cluster"]
    params:
        pos_col       = CFG["options"]["pos_column"],
        dist_method   = CFG["options"]["dist_method"],
        hclust_method = CFG["options"]["hclust_method"],
        h_min         = CFG["options"]["h_min"],
        h_max         = CFG["options"]["h_max"]
    script:
        "src/R/cluster_foci.R"


# Gather per-chromosome foci into a single table per sample_set.
rule _mfR_aggregate:
    input:
        tsv = _mfR_gather_chroms
    output:
        tsv = CFG["dirs"]["foci"] + "{sample_set}/{sample_set}.foci.tsv"
    run:
        import pandas as pd
        frames = []
        for f in input.tsv:
            df = pd.read_csv(f, sep = "\t")
            if df.shape[0] > 0:
                frames.append(df)
        if frames:
            out = pd.concat(frames, ignore_index = True)
            if {"Chromosome", "group"}.issubset(out.columns):
                out["focus_id"] = out["Chromosome"].astype(str) + ":" + out["group"].astype(str)
        else:
            out = pd.read_csv(input.tsv[0], sep = "\t")  # all chromosomes empty -> header only
        out.to_csv(output.tsv, sep = "\t", index = False)


# Symlinks the final output table into the module results directory (under '99-outputs/')
rule _mfR_output_tsv:
    input:
        tsv = str(rules._mfR_aggregate.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{sample_set}.foci.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _mfR_all:
    input:
        expand(
            str(rules._mfR_output_tsv.output.tsv),
            sample_set = CFG["sample_set"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
