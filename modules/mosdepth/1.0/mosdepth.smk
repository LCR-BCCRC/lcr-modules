#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giuliano Banco
# Module Author:    Giuliano Banco
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import sys

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
# `CFG` is a shortcut to `config["lcr-modules"]["mosdepth"]`
CFG = op.setup_module(
    name = "mosdepth",
    version = "1.0",
    subdirectories = ["inputs", "mosdepth", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mosdepth_input_bam,
    _mosdepth_output_common,
    _mosdepth_all,


##### OPTIONS #####

NO_PER_BASE = CFG["options"].get("no_per_base", False)
BY_VALUE = CFG["options"].get("by", None)
MOSDEPTH_ARGS = CFG["options"].get("mosdepth_args", "")

# Avoid letting users hide output-changing flags inside mosdepth_args
for bad_flag in ["--no-per-base", "-n", "--by", "-b"]:
    if bad_flag in MOSDEPTH_ARGS.split():
        sys.exit(
            f"Do not pass {bad_flag} via mosdepth_args. "
            "Use config['lcr-modules']['mosdepth']['options']['no_per_base'] "
            "and ['by'] so declared outputs stay in sync."
        )

PREFIX = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output"

INTERNAL_GLOBAL_DIST = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.mosdepth.global.dist.txt"
INTERNAL_SUMMARY = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.mosdepth.summary.txt"
INTERNAL_PER_BASE_BED = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.per-base.bed.gz"
INTERNAL_PER_BASE_CSI = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.per-base.bed.gz.csi"
INTERNAL_REGIONS_BED = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.regions.bed.gz"
INTERNAL_REGION_DIST = CFG["dirs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/output.mosdepth.region.dist.txt"

MOSDEPTH_OUTPUTS = [
    INTERNAL_GLOBAL_DIST,
    INTERNAL_SUMMARY,
]

if not NO_PER_BASE:
    MOSDEPTH_OUTPUTS.extend([
        INTERNAL_PER_BASE_BED,
        INTERNAL_PER_BASE_CSI,
    ])

if BY_VALUE not in [None, "", False]:
    MOSDEPTH_OUTPUTS.extend([
        INTERNAL_REGIONS_BED,
        INTERNAL_REGION_DIST,
    ])

COMMON_TARGETS = []
OPTIONAL_TARGETS = []


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mosdepth_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Get BAM coverage using Mosdepth
rule _mosdepth:
    input:
        bam = str(rules._mosdepth_input_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        MOSDEPTH_OUTPUTS
    log:
        stdout = CFG["logs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/mosdepth.stdout.log",
        stderr = CFG["logs"]["mosdepth"] + "{seq_type}--{genome_build}/{sample_id}/mosdepth.stderr.log"
    params:
        prefix = PREFIX,
        mosdepth_args = MOSDEPTH_ARGS,
        no_per_base_arg = "--no-per-base" if NO_PER_BASE else "",
        by_arg = f"--by {BY_VALUE}" if BY_VALUE not in [None, "", False] else ""
    conda:
        CFG["conda_envs"]["mosdepth"]
    threads:
        CFG["threads"]["mosdepth"]
    resources:
        **CFG["resources"]["mosdepth"]
    shell:
        op.as_one_line("""
        mosdepth
            {params.mosdepth_args}
            {params.no_per_base_arg}
            {params.by_arg}
            --fasta {input.fasta}
            --threads {threads}
            {params.prefix}
            {input.bam}
            > {log.stdout}
            2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mosdepth_output_common:
    input:
        global_dist = INTERNAL_GLOBAL_DIST,
        summary = INTERNAL_SUMMARY
    output:
        global_dist = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}.output.mosdepth.global.dist.txt",
        summary = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}.output.mosdepth.summary.txt"
    run:
        op.relative_symlink(input.global_dist, output.global_dist, in_module=True)
        op.relative_symlink(input.summary, output.summary, in_module=True)


COMMON_TARGETS.extend([
    str(rules._mosdepth_output_common.output.global_dist),
    str(rules._mosdepth_output_common.output.summary),
])


if not NO_PER_BASE:
    localrules:
        _mosdepth_output_per_base,

    rule _mosdepth_output_per_base:
        input:
            per_base_bed = INTERNAL_PER_BASE_BED,
            per_base_csi = INTERNAL_PER_BASE_CSI
        output:
            per_base_bed = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}.output.per-base.bed.gz",
            per_base_csi = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}.output.per-base.bed.gz.csi"
        run:
            op.relative_symlink(input.per_base_bed, output.per_base_bed, in_module=True)
            op.relative_symlink(input.per_base_csi, output.per_base_csi, in_module=True)

    OPTIONAL_TARGETS.extend([
        str(rules._mosdepth_output_per_base.output.per_base_bed),
        str(rules._mosdepth_output_per_base.output.per_base_csi),
    ])


if BY_VALUE not in [None, "", False]:
    localrules:
        _mosdepth_output_regions,

    rule _mosdepth_output_regions:
        input:
            regions_bed = INTERNAL_REGIONS_BED,
            region_dist = INTERNAL_REGION_DIST
        output:
            regions_bed = CFG["dirs"]["outputs"] + "bed/{seq_type}--{genome_build}/{sample_id}.output.regions.bed.gz",
            region_dist = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{sample_id}.output.mosdepth.region.dist.txt"
        run:
            op.relative_symlink(input.regions_bed, output.regions_bed, in_module=True)
            op.relative_symlink(input.region_dist, output.region_dist, in_module=True)

    OPTIONAL_TARGETS.extend([
        str(rules._mosdepth_output_regions.output.regions_bed),
        str(rules._mosdepth_output_regions.output.region_dist),
    ])


# Generates the target sentinels for each run, which generate the symlinks
rule _mosdepth_all:
    input:
        expand(
            COMMON_TARGETS + OPTIONAL_TARGETS,
            zip,
            seq_type = CFG["samples"]["seq_type"],
            genome_build = CFG["samples"]["genome_build"],
            sample_id = CFG["samples"]["sample_id"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)