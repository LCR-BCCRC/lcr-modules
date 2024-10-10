#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
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
# `CFG` is a shortcut to `config["lcr-modules"]["freebayes"]`
CFG = op.setup_module(
    name = "freebayes",
    version = "1.0",
    subdirectories = ["inputs", "freebayes", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _freebayes_input_bam,
    _freebayes_output_vcf,
    _freebayes_all,

##### CROSSMAP FUNCTIONS #####

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _freebayes_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"], 
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai", 
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    group: 
        "input_and_run"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _freebayes_run:
    input:
        bam = str(rules._freebayes_input_bam.output.bam),
        bai = str(rules._freebayes_input_bam.output.bai), 
        crai = str(rules._freebayes_input_bam.output.crai),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["freebayes"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz", 
        tbi = CFG["dirs"]["freebayes"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["freebayes"] + "{seq_type}--{genome_build}/{sample_id}/freebayes.stderr.log"
    params:
        opts = CFG["options"]["freebayes"]
    conda:
        CFG["conda_envs"]["freebayes"]
    group: 
        "input_and_run"
    threads:
        CFG["threads"]["freebayes"]
    resources:
        **CFG["resources"]["freebayes"]    # All resources necessary can be included and referenced from the config files.
    shell:
        op.as_one_line("""
            freebayes-parallel <(fasta_generate_regions.py {input.fasta}.fai 100000) {threads} -f {input.fasta} {input.bam} | bgzip -o {output.vcf} 2>> {log.stderr} &&
            tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _freebayes_output_vcf:
    input:
        vcf = str(rules._freebayes_run.output.vcf), 
        tbi = str(rules._freebayes_run.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.freebayes.vcf.gz", 
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.freebayes.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _freebayes_all:
    input:
        expand(
            [
                str(rules._freebayes_output_vcf.output.vcf),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
