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
# `CFG` is a shortcut to `config["lcr-modules"]["cutadapt"]`
CFG = op.setup_module(
    name = "cutadapt",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "fastqc_before", "cutadapt", "fastqc_after", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _cutadapt_input_fastq,
    _cutadapt_output_fastq,
    _cutadapt_all,

sample_ids_cutadapt = list(CFG['samples']['sample_id'])

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _cutadapt_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    group: "cutadapt"
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)

rule _cutadapt_fastqc_before:
    input:
        fastq_1 = str(rules._cutadapt_input_fastq.output.fastq_1),
        fastq_2 = str(rules._cutadapt_input_fastq.output.fastq_2)
    output:
        report_1 = CFG["dirs"]["fastqc_before"] + "fastq/{seq_type}/{sample_id}.R1_fastqc.html",
        report_2 = CFG["dirs"]["fastqc_before"] + "fastq/{seq_type}/{sample_id}.R2_fastqc.html"
    log:
        stdout = CFG["logs"]["fastqc_before"] + "{seq_type}/{sample_id}/fastqc_before.stdout.log",
        stderr = CFG["logs"]["fastqc_before"] + "{seq_type}/{sample_id}/fastqc_before.stderr.log"
    params:
        opts = CFG["options"]["fastqc"]
    conda:
        CFG["conda_envs"]["cutadapt"]
    threads:
        CFG["threads"]["fastqc"]
    resources:
        **CFG["resources"]["fastqc"]
    group: "cutadapt"
    wildcard_constraints:
        sample_id = "|".join(sample_ids_cutadapt)
    shell:
        op.as_one_line("""
        fastqc
        {input.fastq_1}
        {input.fastq_2}
        -o $(dirname {output.report_1})/
        {params.opts}
        -t {threads}
        > {log.stdout}
        2> {log.stderr}
        """)

# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _cutadapt_run:
    input:
        fastq_1 = str(rules._cutadapt_input_fastq.output.fastq_1),
        fastq_2 = str(rules._cutadapt_input_fastq.output.fastq_2)
    output:
        trimmed_1 = CFG["dirs"]["cutadapt"] + "fastq/{seq_type}/{sample_id}_trimmed.R1.fastq.gz",
        trimmed_2 = CFG["dirs"]["cutadapt"] + "fastq/{seq_type}/{sample_id}_trimmed.R2.fastq.gz"
    log:
        stdout = CFG["logs"]["cutadapt"] + "{seq_type}/{sample_id}/cutadapt.stdout.log",
        stderr = CFG["logs"]["cutadapt"] + "{seq_type}/{sample_id}/cutadapt.stderr.log"
    params:
        opts = CFG["options"]["cutadapt"],
        forward_a = CFG["options"]["forward_a"],
        reverse_a = CFG["options"]["reverse_a"]
    conda:
        CFG["conda_envs"]["cutadapt"]
    threads:
        CFG["threads"]["cutadapt"]
    resources:
        **CFG["resources"]["cutadapt"]
    group: "cutadapt"
    wildcard_constraints:
        sample_id = "|".join(sample_ids_cutadapt)
    shell:
        op.as_one_line("""
        cutadapt
        -j {threads}
        {params.opts}
        -a {params.forward_a}
        -A {params.reverse_a}
        -o {output.trimmed_1}
        -p {output.trimmed_2}
        {input.fastq_1}
        {input.fastq_2}
        > {log.stdout}
        2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _cutadapt_fastqc_after:
    input:
        fastq_1 = str(rules._cutadapt_run.output.trimmed_1),
        fastq_2 = str(rules._cutadapt_run.output.trimmed_2)
    output:
        report_1 = CFG["dirs"]["fastqc_after"] + "fastq/{seq_type}/{sample_id}_trimmed.R1_fastqc.html",
        report_2 = CFG["dirs"]["fastqc_after"] + "fastq/{seq_type}/{sample_id}_trimmed.R2_fastqc.html"
    log:
        stdout = CFG["logs"]["fastqc_after"] + "{seq_type}/{sample_id}/fastqc_after.stdout.log",
        stderr = CFG["logs"]["fastqc_after"] + "{seq_type}/{sample_id}/fastqc_after.stderr.log"
    params:
        opts = CFG["options"]["fastqc"]
    conda:
        CFG["conda_envs"]["cutadapt"]
    threads:
        CFG["threads"]["fastqc"]
    resources:
        **CFG["resources"]["fastqc"]
    group: "cutadapt"
    wildcard_constraints:
        sample_id = "|".join(sample_ids_cutadapt)
    shell:
        op.as_one_line("""
        fastqc
        {input.fastq_1}
        {input.fastq_2}
        -o $(dirname {output.report_1})/
        {params.opts}
        -t {threads}
        > {log.stdout}
        2> {log.stderr}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _cutadapt_output_fastq:
    input:
        report_before_1 = str(rules._cutadapt_fastqc_before.output.report_1),
        trimmed_1 = str(rules._cutadapt_run.output.trimmed_1),
        report_after_1 = str(rules._cutadapt_fastqc_after.output.report_1),
        report_before_2 = str(rules._cutadapt_fastqc_before.output.report_2),
        trimmed_2 = str(rules._cutadapt_run.output.trimmed_2),
        report_after_2 = str(rules._cutadapt_fastqc_after.output.report_2)
    output:
        report_before_1 = CFG["dirs"]["outputs"] + "report/{seq_type}/{sample_id}.R1_fastqc_before.html",
        trimmed_1 = CFG["dirs"]["outputs"] + "fastq/{seq_type}/{sample_id}_trimmed.R1.fastq.gz",
        report_after_1 = CFG["dirs"]["outputs"] + "report/{seq_type}/{sample_id}.R1_fastqc_after.html",
        report_before_2 = CFG["dirs"]["outputs"] + "report/{seq_type}/{sample_id}.R2_fastqc_before.html",
        trimmed_2 = CFG["dirs"]["outputs"] + "fastq/{seq_type}/{sample_id}_trimmed.R2.fastq.gz",
        report_after_2 = CFG["dirs"]["outputs"] + "report/{seq_type}/{sample_id}.R2_fastqc_after.html"
    run:
        op.relative_symlink(input.report_before_1, output.report_before_1, in_module = True),
        op.relative_symlink(input.trimmed_1, output.trimmed_1, in_module = True),
        op.relative_symlink(input.report_after_1, output.report_after_1, in_module = True),
        op.relative_symlink(input.report_before_2, output.report_before_2, in_module = True),
        op.relative_symlink(input.trimmed_2, output.trimmed_2, in_module = True),
        op.relative_symlink(input.report_after_2, output.report_after_2, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _cutadapt_all:
    input:
        expand(
            [
                str(rules._cutadapt_output_fastq.output.report_before_1),
                str(rules._cutadapt_output_fastq.output.trimmed_1),
                str(rules._cutadapt_output_fastq.output.report_after_1),
                str(rules._cutadapt_output_fastq.output.report_before_2),
                str(rules._cutadapt_output_fastq.output.trimmed_2),
                str(rules._cutadapt_output_fastq.output.report_after_2)
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"]
            )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
