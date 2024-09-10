#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jacky Yiu
# Module Author:    Jacky Yiu
# Contributors:     Sierra Gillis


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
from datetime import datetime
import numpy as np

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
# `CFG` is a shortcut to `config["lcr-modules"]["fishhook"]`
CFG = op.setup_module(
    name = "fishhook",
    version = "1.2",
    subdirectories = ["inputs", "prepare_maf", "fishhook", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _fishhook_input_maf,
    _fishhook_input_subsetting_categories,
    _fishhook_prepare_maf,
    _fishhook_install,
    _fishhook_output_tsv,
    _fishhook_aggregate,
    _fishhook_all,


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS =  os.path.abspath(config["lcr-modules"]["fishhook"]["prepare_mafs"])

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _fishhook_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/{sample_set}--{launch_date}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Symlinks the subsetting categories input file into the module results directory (under '00-inputs/')
rule _fishhook_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

# Prepare the maf file for the input to fishHook
checkpoint _fishhook_prepare_maf:
    input:
        maf = expand(
                    str(rules._fishhook_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._fishhook_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/done"
    log:
        CFG["logs"]["prepare_maf"] + "{sample_set}--{launch_date}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "fishHook",
        metadata_cols = CFG["samples"],
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS


# Install fishHook and required R pacakges
# only available from github, not through conda/CRAN/Biocmanager
rule _fishhook_install:
    output:
        complete = CFG["dirs"]["inputs"] + "fishhook_installed.success"
    conda:
        CFG["conda_envs"]["fishhook"]
    log:
        input = CFG["logs"]["inputs"] + "install_fishhook.log"
    shell:
        """
        R -q --vanilla -e 'devtools::install_github("mskilab/gUtils")' >> {log.input} &&
        R -q --vanilla -e 'devtools::install_github("mskilab/gTrack")' >> {log.input} &&
        R -q --vanilla -e 'devtools::install_github("mskilab/fishHook")' >> {log.input} &&
        touch {output.complete}
        """

# Get gene list input only if that method is specified in the yaml (instead of using tiles)
def get_input_if_gene_mode(wildcards):
    if config["lcr-modules"]["fishhook"]["options"]["target_gene_list"]:
        return reference_files("downloads/gencode-33/gencode.annotation.grch37.gtf")
    else:
        return ""

# Actual fishHook run
rule _fishhook_run:
    input:
        fishhook = ancient(str(rules._fishhook_install.output.complete)),
        maf = CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/{md5sum}.maf",
        content = CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/{md5sum}.maf.content"
    output:
        tsv = CFG["dirs"]["fishhook"] + "{sample_set}--{launch_date}/{md5sum}.fishhook.tsv"
    conda:
        CFG["conda_envs"]["fishhook"]
    log:
        log = CFG["logs"]["fishhook"] + "{sample_set}--{launch_date}--{md5sum}_run_fishook.log"
    threads:
        CFG["threads"]["fishhook"]
    resources:
        **CFG["resources"]["fishhook"]
    params:
        include_silent = CFG["options"]["include_silent_mutation"],
        tiles_size = CFG["options"]["tiles_size"],
        target_gene_list = CFG["options"]["target_gene_list"],
        gene_list = get_input_if_gene_mode,
        target_gene_list_only_protein_coding = CFG["options"]["target_gene_list_only_protein_coding"],
        covariates = CFG["options"]["covariates"]
    script:
        "scr/R/run_fishhook.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _fishhook_output_tsv:
    input:
        tsv = str(rules._fishhook_run.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{sample_set}--{launch_date}/{md5sum}.fishhook.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["fishhook"]
    checkpoint_output = os.path.dirname(str(checkpoints._fishhook_prepare_maf.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum}.maf.content")
    return expand(
        [
            CFG["dirs"]["outputs"] + "tsv/{{sample_set}}--{{launch_date}}/{md5sum}.fishhook.tsv"
        ],
        md5sum = SUMS
        )

# Aggregates outputs to remove md5sum from rule all
rule _fishhook_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{sample_set}--{launch_date}.done"
    shell:
        op.as_one_line("""touch {output.aggregate}""")

# Generates the target sentinels for each run, which generate the symlinks
rule _fishhook_all:
    input:
        expand(
            [
                CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/done",
                str(rules._fishhook_aggregate.output.aggregate),
            ],
            sample_set=CFG["sample_set"],
            launch_date = launch_date)

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
