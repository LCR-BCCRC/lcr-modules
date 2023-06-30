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
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "fishhook", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _fishhook_input_maf,
    _install_fishhook,
    _run_fishhook,
    _fishhook_output_tsv,
    _fishhook_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _fishhook_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

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
        R -q -e 'devtools::install_github(c("jokergoo/ComplexHeatmap","mskilab/gUtils","mskilab/gTrack","mskilab/gChain", "mskilab/skitools", "mskilab/fishHook"))' >> {log.input} &&
        touch {output.complete}"""


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _run_fishhook:
    input:
        fishhook = CFG["dirs"]["inputs"] + "fishhook_installed.success",
        maf = rules._fishhook_input_maf.output.maf,
    output:
        tsv = CFG["dirs"]["fishhook"] + "output.tsv"
    log:
        stdout = CFG["logs"]["fishhook"] + "run_fishhook.stdout.log",
        stderr = CFG["logs"]["fishhook"] + "run_fishhook.stderr.log"
    params:
        tiles_size = CFG["options"]["tiles_size"]
    conda:
        CFG["conda_envs"]["fishhook"]
    threads:
        CFG["threads"]["fishhook"]
    resources:
        **CFG["resources"]["fishhook"]    # All resources necessary can be included and referenced from the config files.
    script:
        "scr/R/run_fishhook.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _fishhook_output_tsv:
    input:
        tsv = str(rules._run_fishhook.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/output.fishhook.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _fishhook_all:
    input:
        tsv = rules._fishhook_output_tsv.output.tsv




##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
