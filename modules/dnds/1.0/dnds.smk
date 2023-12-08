#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     N/A


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
# `CFG` is a shortcut to `config["lcr-modules"]["dnds"]`
CFG = op.setup_module(
    name = "dnds",
    version = "1.0",
    subdirectories = ["inputs", "dnds", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _dnds_input_maf,
    _dnds_input_subsets,
    _dnds_prepare_maf,
    _install_dnds,
    _dnds_output_tsv,
    _dnds_aggregate,
    _dnds_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS =  os.path.abspath(config["lcr-modules"]["dnds"]["prepare_mafs"])

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _dnds_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _dnds_input_subsets:
    input:
        sample_sets = CFG["inputs"]["sample_sets"]
    output:
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)


# Prepare the maf file for the input to MutSig2CV
checkpoint _dnds_prepare_maf:
    input:
        maf = expand(
                    str(rules._dnds_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        sample_sets = str(rules._dnds_input_subsets.output.sample_sets)
    output:
        CFG["dirs"]["inputs"] + "{sample_set}--{launch_date}/done"
    log:
        CFG["logs"]["inputs"] + "{sample_set}--{launch_date}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "dNdS",
        seq_type = CFG["samples"]["seq_type"].unique(),
        metadata = CFG["samples"][["sample_id","seq_type","genome_build","cohort","pathology","unix_group","time_point"]].to_numpy(na_value='')
    script:
        PREPARE_MAFS


# Install dNdS
# only available from github, not through conda/CRAN/Biocmanager
rule _install_dnds:
    output:
        complete = CFG["dirs"]["inputs"] + "dnds_installed.success"
    conda:
        CFG["conda_envs"]["dnds"]
    log:
        input = CFG["logs"]["inputs"] + "install_dnds.log"
    shell:
        """
        R -q -e 'devtools::install_github("im3sanger/dndscv")' >> {log.input} &&
        touch {output.complete}"""


# Actual dNdS run
rule _dnds_run:
    input:
        dnds = ancient(str(CFG["dirs"]["inputs"] + "dnds_installed.success")),
        maf = CFG["dirs"]["inputs"] + "{sample_set}--{launch_date}/{md5sum}.maf",
        content = CFG["dirs"]["inputs"] + "{sample_set}--{launch_date}/{md5sum}.maf.content"
    output:
        dnds_sig_genes = CFG["dirs"]["dnds"] + "{sample_set}--{launch_date}/{md5sum}_sig_genes.tsv",
        annotmuts = CFG["dirs"]["dnds"] + "{sample_set}--{launch_date}/{md5sum}_annotmuts.tsv"
    conda:
        CFG["conda_envs"]["dnds"]
    threads:
        CFG["threads"]["dnds"]
    resources:
        **CFG["resources"]["dnds"]
    params:
        target_genes = CFG["options"]["target_genes"],
        max_muts_per_gene_per_sample = CFG["options"]["max_muts_per_gene_per_sample"]
    script:
        "src/R/run_dnds.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _dnds_output_tsv:
    input:
        tsv = str(rules._dnds_run.output.dnds_sig_genes)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{sample_set}--{launch_date}/{md5sum}.sig_genes.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["dnds"]
    checkpoint_output = os.path.dirname(str(checkpoints._dnds_prepare_maf.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum}.maf.content")
    return expand(
        [
            CFG["dirs"]["outputs"] + "tsv/{{sample_set}}--{{launch_date}}/{md5sum}.sig_genes.tsv"
        ],
        md5sum = SUMS
        )

# aggregates outputs to remove md5sum from rule all
rule _dnds_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{sample_set}--{launch_date}.done"
    shell:
        op.as_one_line("""touch {output.aggregate}""")

# Generates the target sentinels for each run, which generate the symlinks
rule _dnds_all:
    input:
        expand(
            [
                CFG["dirs"]["inputs"] + "{sample_set}--{launch_date}/done",
                str(rules._dnds_aggregate.output.aggregate)
            ],
            sample_set=CFG["sample_set"],
            launch_date = launch_date)

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
