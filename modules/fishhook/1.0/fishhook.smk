#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jacky Yiu
# Module Author:    Jacky Yiu
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
# `CFG` is a shortcut to `config["lcr-modules"]["fishhook"]`
CFG = op.setup_module(
    name = "fishhook",
    version = "1.0",
    subdirectories = ["inputs", "fishhook", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _fishhook_input_maf,
    _fishhook_input_subsets,
    _fishhook_prepare_maf,
    _install_fishhook,
    _run_fishhook,
    _fishhook_output_tsv,
    _fishhook_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _fishhook_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _fishhook_input_subsets:
    input:
        sample_sets = CFG["inputs"]["sample_sets"]
    output:
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)

# Prepare the maf file for the input to Fishhook
rule _fishhook_prepare_maf:
    input:
        maf = expand(
                    str(rules._fishhook_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["seq_types"]
                    ),
        sample_sets = ancient(str(rules._fishhook_input_subsets.output.sample_sets))
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{sample_set}.maf"),
        contents = CFG["dirs"]["inputs"] + "maf/{sample_set}.maf.content"
    log:
        stdout = CFG["logs"]["inputs"] + "{sample_set}/prepare_maf.stdout.log",
        stderr = CFG["logs"]["inputs"] + "{sample_set}/prepare_maf.stderr.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        script = CFG["prepare_mafs"]
    shell:
        op.as_one_line("""
        Rscript {params.script}
        {input.maf}
        {input.sample_sets}
        $(dirname {output.maf})/
        {wildcards.sample_set}
        FishHook
        {params.include_non_coding}
        > {log.stdout} 2> {log.stderr}
        """)


# Install fishhook
# only available from github, not through conda/CRAN/Biocmanager
rule _install_fishhook:
    output:
        complete = CFG["dirs"]["inputs"] + "fishhook_installed.success"
    conda:
        CFG["conda_envs"]["fishhook"]
    log:
        input = CFG["logs"]["inputs"] + "install_fishhook.log"
    shell:
        """
        R -q -e 'devtools::install_github(c("jokergoo/ComplexHeatmap","mskilab/gTrack", "mskilab/fishHook"))' >> {log.input} &&
        touch {output.complete}"""


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _run_fishhook:
    input:
        fishhook = ancient(str(CFG["dirs"]["inputs"] + "fishhook_installed.success")),
        maf = str(rules._fishhook_prepare_maf.output.maf)
    output:
        tsv = CFG["dirs"]["fishhook"] + "{sample_set}/fishhook.output.maf"
    conda:
        CFG["conda_envs"]["fishhook"]
    threads:
        CFG["threads"]["fishhook"]
    resources:
        **CFG["resources"]["fishhook"]     
    params:
        tiles_size = CFG["options"]["tiles_size"],
        coveriate = CFG["options"]["coveriates"],
        include_silent = CFG["options"]["include_silent_mutation"],
        gene_list = CFG["options"]["gene_list"],
        gene_list_only_protein_coding = CFG["options"]["gene_list_only_protein_coding"]
    script:
        "scr/R/run_fishhook.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _fishhook_output_tsv:
    input:
        tsv = str(rules._run_fishhook.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{sample_set}/{sample_set}.fishhook.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _fishhook_all:
    input:
        expand(
            [
                str(rules._fishhook_output_tsv.output.tsv),
            ],
            sample_set=CFG["sample_set"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
