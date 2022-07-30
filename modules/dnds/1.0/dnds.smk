#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
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
# `CFG` is a shortcut to `config["lcr-modules"]["dnds"]`
CFG = op.setup_module(
    name = "dnds",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "dnds", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _dnds_input_maf,
    _dnds_step_2,
    _dnds_output_tsv,
    _dnds_all,


##### RULES #####


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
        sample_sets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_sets, output.sample_sets)


# Prepare the maf file for the input to MutSig2CV
rule _dnds_prepare_maf:
    input:
        maf = expand(
                    str(rules._dnds_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["seq_types"]
                    ),
        sample_sets = ancient(str(rules._dnds_input_subsets.output.sample_sets))
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
        dNdS
        {params.include_non_coding}
        > {log.stdout} 2> {log.stderr}
        """)


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
        maf = str(rules._dnds_prepare_maf.output.maf)
    output:
        dnds_sig_genes = CFG["dirs"]["dnds"] + "{sample_set}/sig_genes.tsv",
        annotmuts = CFG["dirs"]["dnds"] + "{sample_set}/annotmuts.tsv"
    conda:
        CFG["conda_envs"]["dnds"]
    threads:
        CFG["threads"]["dnds"]
    resources:
        **CFG["resources"]["dnds"]
    params:
        running_directory = CFG["dirs"]["mcr"]
    script:
        "src/R/run_dnds.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _dnds_output_tsv:
    input:
        tsv = str(rules._dnds_run.output.dnds_sig_genes)
    output:
        txt = CFG["dirs"]["outputs"] + "tsv/{sample_set}/{sample_set}.sig_genes.tsv"
    run:
        op.relative_symlink(input.txt, output.txt, in_module= True)


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _dnds_output_tsv:
    input:
        tsv = str(rules._dnds_step_2.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _dnds_all:
    input:
        expand(
            [
                str(rules._dnds_output_tsv.output.tsv),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
