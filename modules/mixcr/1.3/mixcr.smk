#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Anita dos Santos
# Module Author:    Manuela Cruz
# Contributors:     Laura Hilton


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "mixcr",
    version = "1.3",
    subdirectories = ["inputs", "mixcr", "outputs"],
)

RECEPTORS = CFG["receptors"]

ig_type = ["IGH","IGK","IGL"]
tr_type = ["TRA", "TRB", "TRD", "TRG"]

RANGE = ig_type + tr_type

if isinstance(RECEPTORS, str):
    RECEPTORS = RECEPTORS.split(" ")
if isinstance(RECEPTORS, list) and len(RECEPTORS) == 1:
    if RECEPTORS[0] == "ALL":
        RECEPTORS = ig_type + tr_type
    if RECEPTORS[0] == "BCR":
        RECEPTORS = ig_type
    if RECEPTORS[0] == "TCR":
        RECEPTORS = tr_type

assert all(receptor in RANGE for receptor in RECEPTORS), (
    "Config 'receptors' value is necessary for specifying target receptors. "
    "Ensure desired receptors are included in config and uppercase. "
    "Choose from: 'ALL', 'BCR', 'TCR' or list of IGH, IGK, IGL, TRA, TRB, TRD, TRG. "
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mixcr_input_fastq,
    _mixcr_output_txt,
    _mixcr_to_fasta,
    _mixcr_output_fasta,
    _mixcr_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mixcr_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)

# Run MiXCR
rule _mixcr_run:
    input:
        fastq_1 = str(rules._mixcr_input_fastq.output.fastq_1),
        fastq_2 = str(rules._mixcr_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"], # Prevent premature deletion of fastqs marked as temp
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
    output:
        results = expand(CFG["dirs"]["mixcr"] + "{{seq_type}}/{{sample_id}}/mixcr.{{sample_id}}.clones_{chain}.tsv", chain = RECEPTORS)
    log:
        stdout = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_run.stdout.log",
        stderr = CFG["logs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr_run.stderr.log"
    resources:
        **CFG["resources"]["mixcr_run"]
    params:
        preset       = op.switch_on_wildcard("seq_type", CFG["options"]["preset"]),
        opts         = op.switch_on_wildcard("seq_type", CFG["options"]["mixcr_run"]),
        prefix       = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}",
        license_file = CFG["inputs"]["license_file"],
        jvmheap      = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        chains       = " ".join(RECEPTORS)
    conda: CFG["conda_envs"]["java"]
    container:
        CFG["container_envs"]["java"]
    threads:
        CFG["threads"]["mixcr_run"]
    shell:
        op.as_one_line("""
        MI_LICENSE_FILE={params.license_file}
        mixcr -Xmx{params.jvmheap}m analyze {params.preset}
        --species hsa
        --assemble-longest-contigs
        -t {threads} {params.opts}
        {input.fastq_1} {input.fastq_2}
        {params.prefix} > {log.stdout} 2> {log.stderr}
        &&
        for chain in {params.chains}; do
            touch {params.prefix}.clones_$chain.tsv ;
        done
        """)

# Convert MiXCR clonotype table to FASTA format for downstream use (e.g. igblast module)
rule _mixcr_to_fasta:
    input:
        mixcr_results = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clones_{chain}.tsv",
        script = CFG["scripts"]["mixcr2fasta"]
    output:
        fasta = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clones_{chain}.VDJseq.fasta",
        seq_info = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clones_{chain}.regions.txt"
    shell:
        op.as_one_line("""
        if [ -s {input.mixcr_results} ]; then
            {input.script} -i {input.mixcr_results} -o {output.fasta} -s {output.seq_info} ;
        else
            touch {output.fasta} {output.seq_info} ;
        fi
        """)

# Symlinks the final clonotype txt files into the module results directory (under '99-outputs/')
rule _mixcr_output_txt:
    input:
        results = CFG["dirs"]["mixcr"] + "{seq_type}/{sample_id}/mixcr.{sample_id}.clones_{chain}.tsv"
    output:
        results = CFG["dirs"]["outputs"] + "txt/{seq_type}/mixcr.{sample_id}.clones_{chain}.tsv"
    wildcard_constraints:
        chain = '[A-Z]+'
    run:
        op.relative_symlink(input.results, output.results, in_module=True)

# Symlinks the per-chain FASTA and regions.txt into '99-outputs/' for use by downstream modules
rule _mixcr_output_fasta:
    input:
        fasta = str(rules._mixcr_to_fasta.output.fasta),
        seq_info = str(rules._mixcr_to_fasta.output.seq_info)
    output:
        fasta = CFG["dirs"]["outputs"] + "fasta/{seq_type}/mixcr.{sample_id}.clones_{chain}.VDJseq.fasta",
        seq_info = CFG["dirs"]["outputs"] + "seq_info/{seq_type}/mixcr.{sample_id}.clones_{chain}.regions.tsv"
    wildcard_constraints:
        chain = '[A-Z]+'
    run:
        op.relative_symlink(input.fasta, output.fasta, in_module=True)
        op.relative_symlink(input.seq_info, output.seq_info, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _mixcr_all:
    input:
        expand(
            expand(
                [
                    str(rules._mixcr_output_txt.output.results),
                    str(rules._mixcr_output_fasta.output.fasta),
                ],
                zip,
                seq_type=CFG["samples"]["seq_type"],
                sample_id=CFG["samples"]["sample_id"],
                allow_missing=True),
            chain=RECEPTORS)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
