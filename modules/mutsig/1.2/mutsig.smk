#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     Sierra Gillis


##### SETUP #####

# Import package with useful functions for developing analysis modules
import sys, os
from os.path import join
import glob
import hashlib
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
# `CFG` is a shortcut to `config["lcr-modules"]["mutsig"]`
CFG = op.setup_module(
    name = "mutsig",
    version = "1.1",
    subdirectories = ["inputs", "prepare_maf","mcr", "mutsig", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutsig_input_maf,
    _mutsig_input_subsetting_categories,
    _mutsig_prepare_maf,
    _mutsig_download_mutsig,
    _mutsig_download_mcr,
    _mutsig_configure_mcr,
    _mutsig_output_txt,
    _mutsig_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS =  os.path.abspath(config["lcr-modules"]["mutsig"]["prepare_mafs"])

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutsig_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/{sample_set}--{launch_date}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Symlinks the subsetting categories input file into the module results directory (under '00-inputs/')
rule _mutsig_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)


# Prepare the maf file for the input to MutSig2CV
checkpoint _mutsig_prepare_maf:
    input:
        maf = expand(
                    str(rules._mutsig_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._mutsig_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/done"
    log:
        CFG["logs"]["prepare_maf"] + "{sample_set}--{launch_date}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "MutSig2CV",
        metadata_cols = CFG["samples"],
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS


# Download the MutSig2CV
rule _mutsig_download_mutsig:
    output:
        mutsig_compressed = temp(CFG["dirs"]["inputs"] + "MutSig2CV/MutSig2CV.tar.gz"),
        mutsig = CFG["dirs"]["mcr"] + "MutSig2CV/MutSig2CV.success"
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
        wget
        -O {output.mutsig_compressed}
        http://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/MutSig2CV.tar.gz
            &&
        tar -xvf {output.mutsig_compressed} -C $(dirname {output.mutsig})
            &&
        touch {output.mutsig}
        """)

# Download the MCR installer
rule _mutsig_download_mcr:
    output:
        mcr_installer = temp(CFG["dirs"]["mcr"] + "MCR_R2013a_glnxa64_installer.zip"),
        local_mcr = directory(CFG["dirs"]["mcr"] + "local_mcr")
    log:
        stdout = CFG["logs"]["mcr"] + "download_mcr.stdout.log",
        stderr = CFG["logs"]["mcr"] + "download_mcr.stderr.log"
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
        wget
        -O {output.mcr_installer}
        https://ssd.mathworks.com/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip
            &&
        unzip {output.mcr_installer} -d $(dirname {output.mcr_installer})
        > {log.stdout} 2> {log.stderr}
            &&
        mkdir {output.local_mcr}
        """)


# Install local MCR
rule _mutsig_install_mcr:
    input:
        mcr = str(rules._mutsig_download_mcr.output.mcr_installer)
    output:
        local_mcr = CFG["dirs"]["mcr"] + "install.success"
    log:
        stdout = CFG["logs"]["mcr"] + "install_mcr.stdout.log",
        stderr = CFG["logs"]["mcr"] + "install_mcr.stderr.log"
    conda:
        CFG["conda_envs"]["matlab"]
    shell:
        op.as_one_line("""
        $(dirname {input.mcr})/install
        -mode silent
        -agreeToLicense yes
        -destinationFolder $(dirname $(readlink -f {output.local_mcr}))/local_mcr
        > {log.stdout} 2> {log.stderr}
            &&
        touch {output.local_mcr}
        """)

# Obtain the path to the matlab conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open(CFG["conda_envs"]["matlab"], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
print(conda_prefix + "/" + h)
MATLAB = str(conda_prefix + "/" + h)

# Configure local MCR
rule _mutsig_configure_mcr:
    input:
        mcr = str(rules._mutsig_download_mcr.output.local_mcr),
        mcr_installed = str(rules._mutsig_install_mcr.output.local_mcr),
        mutsig = str(rules._mutsig_download_mutsig.output.mutsig)
    output:
        local_mcr = CFG["dirs"]["mcr"] + "configure.success",
        configured = MATLAB + "/lib/configure.success"
    conda:
        CFG["conda_envs"]["matlab"]
    params:
        running_directory = os.path.abspath(CFG["dirs"]["mcr"]),
        scripts_directory = os.path.abspath(CFG["src_dir"]),
        conda_prefix = lambda w: workflow.conda_prefix if workflow.conda_prefix else os.path.abspath(".snakemake/conda")
    shell:
        op.as_one_line("""
        cd {input.mcr}/v81/runtime/glnxa64/
            &&
        ln -sf libmwmclmcrrt.so.8.1 libmwmclmcrrt.so.9.0.1
            &&
        cd ../../bin/glnxa64/
            &&
        ln -sf libboost_system.so.1.49.0 libboost_system.so.1.56.0
            &&
        cd {params.running_directory}
            &&
        touch {output.configured}
            &&
        ln -sf {params.scripts_directory}/bash/run_MutSig2CV.sh .
            &&
        ln -sf $(dirname {output.configured})/libncurses.so.6 libncurses.so.5
            &&
        ln -sf $(dirname {output.configured})/libtinfow.so.6 libtinfow.so.6
            &&
        ln -sf {params.running_directory}/sys/java/jre/glnxa64/jre/lib/amd64/server/libjvm.so libjvm.so
            &&
        ln -sf MutSig2CV/mutsig2cv/reference .
            &&
        touch configure.success
        """)


# Actual MutSig2CV run
rule _mutsig_run:
    input:
        mcr = ancient(str(rules._mutsig_download_mcr.output.local_mcr)),
        mcr_installed = ancient(str(rules._mutsig_install_mcr.output.local_mcr)),
        mcr_configured = ancient(str(rules._mutsig_configure_mcr.output.local_mcr)),
        mutsig = ancient(str(rules._mutsig_download_mutsig.output.mutsig)),
        maf = CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/{md5sum}.maf",
        content = CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/{md5sum}.maf.content"
    output:
        mutsig_maf = temp(CFG["dirs"]["mutsig"] + "{sample_set}--{launch_date}--{md5sum}/final_analysis_set.maf"),
        mutsig_sig_genes = CFG["dirs"]["mutsig"] + "{sample_set}--{launch_date}--{md5sum}/sig_genes.txt",
        success = CFG["dirs"]["mutsig"] + "{sample_set}--{launch_date}--{md5sum}/mutsig.success"
    conda:
        CFG["conda_envs"]["matlab"]
    threads:
        CFG["threads"]["mutsig"]
    resources:
        **CFG["resources"]["mutsig"]
    params:
        running_directory = CFG["dirs"]["mcr"]
    shell:
        op.as_one_line("""
        cd {params.running_directory}
            &&
        ./run_MutSig2CV.sh
        local_mcr/v81
        ../01-prepare_maf/{wildcards.sample_set}--{wildcards.launch_date}/{wildcards.md5sum}.maf
        ../03-mutsig/{wildcards.sample_set}--{wildcards.launch_date}--{wildcards.md5sum}/
        > ../03-mutsig/{wildcards.sample_set}--{wildcards.launch_date}--{wildcards.md5sum}/log
            &&
        touch ../03-mutsig/{wildcards.sample_set}--{wildcards.launch_date}--{wildcards.md5sum}/mutsig.success
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutsig_output_txt:
    input:
        txt = str(rules._mutsig_run.output.mutsig_sig_genes)
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{sample_set}--{launch_date}/{md5sum}.sig_genes.txt"
    run:
        op.relative_symlink(input.txt, output.txt, in_module= True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["mutsig"]
    checkpoint_output = os.path.dirname(str(checkpoints._mutsig_prepare_maf.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum}.maf.content")
    return expand(
        [
            CFG["dirs"]["outputs"] + "txt/{{sample_set}}--{{launch_date}}/{md5sum}.sig_genes.txt"
        ],
        md5sum = SUMS
        )

# Aggregates outputs to remove md5sum from rule all
rule _mutsig_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{sample_set}--{launch_date}.done"
    shell:
        op.as_one_line("""touch {output.aggregate}""")


# Generates the target sentinels for each run, which generate the symlinks
rule _mutsig_all:
    input:
        expand(
            [
                CFG["dirs"]["prepare_maf"] + "{sample_set}--{launch_date}/done",
                str(rules._mutsig_aggregate.output.aggregate),
            ],
            sample_set=CFG["sample_set"],
            launch_date = launch_date)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
