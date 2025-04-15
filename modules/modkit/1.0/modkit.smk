#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Haya Shaalan
# Module Author:    Haya Shaalan
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
# `CFG` is a shortcut to `config["lcr-modules"]["modkit"]`
CFG = op.setup_module(
    name = "modkit",
    version = "1.0",
    subdirectories = ["inputs", "modkit", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _modkit_input_bam,
    _modkit_output_tsv,
    _modkit_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _modkit_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"],
        sample_bai = CFG["inputs"]["sample_bai"]
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        sample_crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.sample_bam, output.sample_bam)
        op.absolute_symlink(input.sample_bai, output.sample_bai)
        op.absolute_symlink(input.sample_bai, output.sample_crai)


rule _modkit_pileup:
    input:
       bam = str(rules._modkit_input_bam.output.sample_bam),
       ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        meth = temp(CFG["dirs"]["modkit"] + "{seq_type}--{genome_build}/{sample_id}_methylation.tsv"),
        meth_gz = CFG["dirs"]["modkit"] + "{seq_type}--{genome_build}/{sample_id}_methylation.tsv.gz"
    log:
        stdout = CFG["logs"]["modkit"] + "{seq_type}--{genome_build}/{sample_id}/pileup.stdout.log",
        stderr = CFG["logs"]["modkit"] + "{seq_type}--{genome_build}/{sample_id}/pileup.stderr.log"
    conda: 
        CFG["conda_envs"]["modkit"]
    params:
        CFG["options"]["header"]
    threads:
        CFG["threads"]["modkit"]
    resources:
        **CFG["resources"]["modkit"]    # All resources necessary can be included and referenced from the config files.
    shell:
        """
        # Run modkit pileup and generate output
        modkit pileup {input.bam} {output.meth} --cpg -t {threads} --ref {input.ref} > {log.stdout} 2> {log.stderr}
        
        # Add header and compress the final output
        (echo -e "chrom\tstart\tend\tmodified_base_code\tscore\tstrand\tstart_position\tend_position\tcolor\tNvalid_cov\tpercent_modified\tNmod\tNcanonical\tNother_mod\tNdelete\tNfail\tNdiff\tNnocall"; cat {output.meth}) | gzip > {output.meth_gz}
        """

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _modkit_output_tsv:
    input:
        tsv_gz = str(rules._modkit_pileup.output.meth_gz)
    output:
        tsv_gz = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{sample_id}.output.tsv.gz"
    run:
        op.relative_symlink(input.tsv_gz, output.tsv_gz, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _modkit_all:
    input:
        expand(
            [
                str(rules._modkit_output_tsv.output.tsv_gz),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
