#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
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
# `CFG` is a shortcut to `config["lcr-modules"]["jabba"]`
CFG = op.setup_module(
    name = "jabba",
    version = "1.0",
    subdirectories = ["inputs", "fragcounter", "outputs"],
)


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _jabba_input_bam,
    _jabba_output_tumour_rds,
    _jabba_output_normal_rds,
    _jabba_all


##### RULES #####


rule _jabba_install_fragcounter:
    output:
        complete = CFG["dirs"]["fragcounter"] + "fragcounter.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        """
        Rscript -e 'remotes::install_github("mskilab/fragCounter")' &&
        touch {output.complete}
        """

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _jabba_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"]+ "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)


rule _jabba_run_fragcounter:
    input:
        installed = str(rules._jabba_install_fragcounter.output.complete),
        bam = str(rules._jabba_input_bam.output.bam),
        gc = reference_files("genomes/{genome_build}/annotations/jabba/gc1000.rds"),
        map = reference_files("genomes/{genome_build}/annotations/jabba/map1000.rds")
    output:
        rds = CFG["dirs"]["fragcounter"] + "{seq_type}--{genome_build}/{sample_id}/cov.rds"
    log:
        stdout = CFG["logs"]["fragcounter"] + "{seq_type}--{genome_build}/{sample_id}/fc.stdout.log",
        stderr = CFG["logs"]["fragcounter"] + "{seq_type}--{genome_build}/{sample_id}/fc.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line(""" 
        export PATH=${{PATH}}:$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))') &&
        frag -b {input.bam} -d `dirname {input.gc}` -w 1000 -o `dirname {output.rds}` > {log.stdout} 2> {log.stderr}
        """)


rule _jabba_output_tumour_rds:
    input:
        rds = CFG["dirs"]["fragcounter"] + "{seq_type}--{genome_build}/{tumour_id}/cov.rds"
    output:
        rds = CFG["dirs"]["outputs"] + "rds/{seq_type}--{genome_build}/tumour/{tumour_id}.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

rule _jabba_output_normal_rds:
    input:
        rds = CFG["dirs"]["fragcounter"] + "{seq_type}--{genome_build}/{normal_id}/cov.rds"
    output:
        rds = CFG["dirs"]["outputs"] + "rds/{seq_type}--{genome_build}/normal/{normal_id}.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Generates the target sentinels for each run, which generate the symlinks
rule _jabba_all:
    input:
        expand(
            [
                str(rules._jabba_output_tumour_rds.output.rds),
                str(rules._jabba_output_normal_rds.output.rds)
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
