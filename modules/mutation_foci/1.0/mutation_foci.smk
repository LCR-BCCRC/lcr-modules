#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giuliano Banco
# Module Author:    Giuliano Banco
# Contributors:     Houman Layegh Mirhosseini


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
# `CFG` is a shortcut to `config["lcr-modules"]["mutation_foci"]`
CFG = op.setup_module(
    name = "mutation_foci",
    version = "1.0",
    subdirectories = ["inputs", "prepare", "foci", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutation_foci_input_maf,
    _mutation_foci_input_subsets,
    _mutation_foci_prepare_maf,
    _mutation_foci_aggregate,
    _mutation_foci_output_tsv,
    _mutation_foci_all,


##### FUNCTIONS #####


def _mutation_foci_gather_chroms(wildcards):
    """Gather the per-chromosome foci tables for one sample_set."""
    return expand(
        CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.foci.tsv",
        sample_set = wildcards.sample_set,
        chrom = CFG["chromosomes"]
    )


##### RULES #####


# Symlinks the input MAF into the module results directory (under '00-inputs/')
rule _mutation_foci_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Symlinks the sample-set membership table into the module results directory
rule _mutation_foci_input_subsets:
    input:
        sample_sets = CFG["inputs"]["sample_sets"]
    output:
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)


# Subset the master MAF to one sample_set and keep only non-coding mutations.
# seq_types listed in the config are combined into a single per-sample_set MAF.
rule _mutation_foci_prepare_maf:
    input:
        maf = expand(
            str(rules._mutation_foci_input_maf.output.maf),
            allow_missing = True,
            seq_type = CFG["seq_types"]
        ),
        sample_sets = ancient(str(rules._mutation_foci_input_subsets.output.sample_sets))
    output:
        maf = temp(CFG["dirs"]["prepare"] + "{sample_set}.noncoding.maf")
    log:
        log = CFG["logs"]["prepare"] + "{sample_set}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["mutation_foci"]
    params:
        coding_classes = CFG["options"]["coding_variant_classifications"],
        sample_id_col = CFG["options"]["sample_id_column"],
        sample_set_col = CFG["options"]["sample_set_column"],
        maf_sample_col = CFG["options"]["maf_sample_column"]
    script:
        "src/R/prepare_maf.R"


# Cluster non-coding mutation positions into "foci", scattered by chromosome.
# This is the heavy step (one job per sample_set x chromosome).
rule _mutation_foci_cluster:
    input:
        maf = str(rules._mutation_foci_prepare_maf.output.maf)
    output:
        tsv  = CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.foci.tsv",
        plot = CFG["dirs"]["foci"] + "{sample_set}/chromosomes/{chrom}.silhouette.pdf"
    log:
        log = CFG["logs"]["foci"] + "{sample_set}/{chrom}.cluster.log"
    conda:
        CFG["conda_envs"]["mutation_foci"]
    threads:
        CFG["threads"]["cluster"]
    resources:
        # dist() is O(n^2) in unique positions; bump mem for dense chromosomes.
        **CFG["resources"]["cluster"]
    params:
        chrom_col = CFG["options"]["chrom_column"],
        pos_col = CFG["options"]["pos_column"],
        dist_method = CFG["options"]["dist_method"],
        hclust_method = CFG["options"]["hclust_method"],
        h_min = CFG["options"]["h_min"],
        h_max = CFG["options"]["h_max"]
    script:
        "src/R/cluster_foci.R"


# Gather per-chromosome foci into a single table per sample_set.
rule _mutation_foci_aggregate:
    input:
        tsv = _mutation_foci_gather_chroms
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
rule _mutation_foci_output_tsv:
    input:
        tsv = str(rules._mutation_foci_aggregate.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{sample_set}.foci.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutation_foci_all:
    input:
        expand(
            str(rules._mutation_foci_output_tsv.output.tsv),
            sample_set = CFG["sample_set"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
