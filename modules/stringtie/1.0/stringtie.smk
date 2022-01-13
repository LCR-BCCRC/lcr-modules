#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Krysta Coyle
# Module Author:    Krysta Coyle
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
# `CFG` is a shortcut to `config["lcr-modules"]["stringtie"]`
CFG = op.setup_module(
    name = "stringtie",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "stringtie", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _stringtie_input_bam,
    _stringtie_output_gtf,
    _stringtie_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _stringtie_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    group: 
        "input_and_step_1"
    run:
        op.absolute_symlink(input.bam, output.bam),
        op.absolute_symlink(input.bam + ".bai", output.bai)


# Run stringtie
rule _stringtie_run:
    input:
        bam = str(rules._stringtie_input_bam.output.bam),
        ref_gtf = reference_files("genomes/{genome_build}/annotations/gencode_annotation-33.gtf"),
        XS_script = CFG["inputs"]["XS_script"]
    output:
        gtf = CFG["dirs"]["stringtie"] + "{seq_type}--{genome_build}/{sample_id}/output.gtf"
    log:
        stdout = CFG["logs"]["stringtie"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["stringtie"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["stringtie_run"]
    conda:
        CFG["conda_envs"]["stringtie"]
    threads:
        CFG["threads"]["stringtie_run"]
    resources:
        **CFG["resources"]["stringtie_run"]    # All resources necessary can be included and referenced from the config files.
    shell:
        # stringtie –l label –p 16 -G hg19gene.gtf –o output.gtf input.bam
        op.as_one_line("""
        samtools view -h {input.bam} | \
        awk -v strType=2 -f {input.XS_script} | \
        stringtie -o {output.gtf} -G {input.ref_gtf} {params.opts} \
        -p {threads} > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _stringtie_output_gtf:
    input:
        gtf = str(rules._stringtie_run.output.gtf)
    output:
        gtf = CFG["dirs"]["outputs"] + "gtf/{seq_type}--{genome_build}/{sample_id}.output.filt.gtf"
    run:
        op.relative_symlink(input.gtf, output.gtf, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _stringtie_all:
    input:
        expand(
            [
                str(rules._stringtie_output_gtf.output.gtf),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
