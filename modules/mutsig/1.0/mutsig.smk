#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import sys, os
from os.path import join
import glob
import hashlib
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
# `CFG` is a shortcut to `config["lcr-modules"]["mutsig"]`
CFG = op.setup_module(
    name = "mutsig",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "mcr", "mutsig", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _mutsig_input_maf,
    _mutsig_step_2,
    _mutsig_output_txt,
    _mutsig_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutsig_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutsig_input_subsets:
    input:
        sample_sets = CFG["inputs"]["sample_sets"]
    output:
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)


# Prepare the maf file for the input to MutSig2CV
rule _mutsig_prepare_maf:
    input:
        maf = expand(
                    str(rules._mutsig_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["seq_types"]
                    ),
        sample_sets = str(rules._mutsig_input_subsets.output.sample_sets)
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
        MutSig2CV
        {params.include_non_coding}
        > {log.stdout} 2> {log.stderr}
        """)


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
        mcr = str(rules._mutsig_download_mcr.output.local_mcr),
        mcr_installed = str(rules._mutsig_install_mcr.output.local_mcr),
        mcr_configured = str(rules._mutsig_configure_mcr.output.local_mcr),
        mutsig = str(rules._mutsig_download_mutsig.output.mutsig),
        maf = str(rules._mutsig_prepare_maf.output.maf)
    output:
        mutsig_maf = temp(CFG["dirs"]["mutsig"] + "{sample_set}/final_analysis_set.maf"),
        mutsig_sig_genes = CFG["dirs"]["mutsig"] + "{sample_set}/sig_genes.txt",
        success = CFG["dirs"]["mutsig"] + "{sample_set}/mutsig.success"
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
        ../00-inputs/maf/{wildcards.sample_set}.maf
        ../02-mutsig/{wildcards.sample_set}/
        > ../02-mutsig/{wildcards.sample_set}/log
            &&
        touch ../02-mutsig/{wildcards.sample_set}/mutsig.success
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _mutsig_output_txt:
    input:
        txt = str(rules._mutsig_run.output.mutsig_sig_genes)
    output:
        txt = CFG["dirs"]["outputs"] + "txt/{sample_set}/{sample_set}.sig_genes.txt"
    run:
        op.relative_symlink(input.txt, output.txt, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutsig_all:
    input:
        expand(
            [
                str(rules._mutsig_output_txt.output.txt),
                # TODO: If applicable, add other output rules here
            ],
            sample_set=CFG["sample_set"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
