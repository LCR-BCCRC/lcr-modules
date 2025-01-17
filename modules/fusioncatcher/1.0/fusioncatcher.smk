#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Daniel Nicorici
# Module Author:    Ryan Morin
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
# `CFG` is a shortcut to `config["lcr-modules"]["fusioncatcher"]`
CFG = op.setup_module(
    name = "fusioncatcher",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "fusioncatcher", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _fusioncatcher_input_bam,
    _fusioncatcher_step_2,
    _fusioncatcher_output_tsv,
    _fusioncatcher_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _fusioncatcher_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    group: 
        "input_and_step_1"
    run:
        op.absolute_symlink(input.bam, output.bam)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _fusioncatcher_step_1:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        tsv = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.tsv"
    log:
        stdout = CFG["logs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        **CFG["resources"]["step_1"]
    group: 
        "input_and_step_1"
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --tumour {input.tumour_bam} --normal {input.normal_bam}
        --ref-fasta {input.fasta} --output {output.tsv} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _fusioncatcher_step_2:
    input:
        tsv = str(rules._fusioncatcher_step_1.output.tsv)
    output:
        tsv = CFG["dirs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.tsv"
    log:
        stderr = CFG["logs"]["fusioncatcher"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.tsv} > {output.tsv} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _fusioncatcher_output_tsv:
    input:
        tsv = str(rules._fusioncatcher_step_2.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _fusioncatcher_all:
    input:
        expand(
            [
                str(rules._fusioncatcher_output_tsv.output.tsv),
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
