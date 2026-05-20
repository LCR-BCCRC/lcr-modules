#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


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
# `CFG` is a shortcut to `config["lcr-modules"]["igseqr"]`
CFG = op.setup_module(
    name = "igseqr",
    version = "1.0",
    subdirectories = ["inputs", "igseqr", "outputs"],
)

assert len(CFG["inputs"]["hisat_ref_url"]) == 1, (
    "Config 'hisat_ref_url' must contain exactly one entry (version_string: url). "
    "To use a different index, replace the key and URL."
)

HISAT_REF_VERSION = list(CFG["inputs"]["hisat_ref_url"].keys())[0]
HISAT_REF_URL     = CFG["inputs"]["hisat_ref_url"][HISAT_REF_VERSION]
HISAT_REF_DIR     = CFG["dirs"]["inputs"] + "hisat_ref"

CHAINS = CFG["options"]["chains"]
if isinstance(CHAINS, str):
    CHAINS = CHAINS.split()

assert all(c in {"IGH", "IGKL"} for c in CHAINS), (
    "Config 'options.chains' must contain 'IGH', 'IGKL', or both."
)

# Map the configured chain list to the igseqr --chain flag.
# --chain B runs both heavy and light chains in one invocation.
if set(CHAINS) >= {"IGH", "IGKL"}:
    CHAIN_ARG = "B"
elif "IGH" in CHAINS:
    CHAIN_ARG = "H"
else:
    CHAIN_ARG = "L"

# Pre-expand chain into output path templates; {seq_type} and {sample_id}
# remain as Snakemake wildcards.
_IGSEQR_RUN_OUTPUTS = expand(
    [
        CFG["dirs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_transcripts.fasta",
        CFG["dirs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_report.tsv",
        CFG["dirs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_dominant_report.tsv",
        CFG["dirs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/{sample_id}_{chain}_TPM_filtered.fasta",
    ],
    chain=CHAINS,
    allow_missing=True,
)

_IGSEQR_OUTPUT_FILES = expand(
    [
        CFG["dirs"]["outputs"] + "fasta/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_transcripts.fasta",
        CFG["dirs"]["outputs"] + "tsv/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_report.tsv",
        CFG["dirs"]["outputs"] + "tsv/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_dominant_report.tsv",
        CFG["dirs"]["outputs"] + "fasta/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_TPM_filtered.fasta",
    ],
    chain=CHAINS,
    allow_missing=True,
)

# Define rules to be run locally when using a compute cluster
localrules:
    _igseqr_input_fastq,
    _igseqr_output_files,
    _igseqr_all,


##### RULES #####


# Downloads the HISAT2 splice-aware index specified in hisat_ref_url.
# The index is stored under {inputs_dir}/hisat_ref/{version}/ and the prefix
# {inputs_dir}/hisat_ref/{version}/genome is passed to igseqr.
rule _igseqr_get_hisat_ref:
    params:
        ref_dir = HISAT_REF_DIR,
        version = HISAT_REF_VERSION,
        url     = HISAT_REF_URL,
    output:
        complete = HISAT_REF_DIR + "/hisat2_" + HISAT_REF_VERSION + "_downloaded.success"
    log:
        stdout = CFG["logs"]["igseqr"] + "get_hisat_ref.stdout.log",
        stderr = CFG["logs"]["igseqr"] + "get_hisat_ref.stderr.log",
    resources:
        **CFG["resources"]["get_hisat_ref"]
    shell:
        '''
        wget -O {params.ref_dir}/{params.version}.tar.gz "{params.url}" \
            > {log.stdout} 2> {log.stderr}
        tar -xzf {params.ref_dir}/{params.version}.tar.gz -C {params.ref_dir}/ \
            >> {log.stdout} 2>> {log.stderr}
        rm {params.ref_dir}/{params.version}.tar.gz
        touch {output.complete}
        '''


# Symlinks the input FASTQ pair into the module results directory (under '00-inputs/')
rule _igseqr_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"],
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R1.fastq.gz",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}/{sample_id}.R2.fastq.gz",
    run:
        op.absolute_symlink(input.fastq_1, output.fastq_1)
        op.absolute_symlink(input.fastq_2, output.fastq_2)


# Run IgSeqR for the configured chain(s) in one invocation.
# CHAIN_ARG is "H", "L", or "B" (both), derived from options.chains at load time.
rule _igseqr_run:
    input:
        fastq_1 = str(rules._igseqr_input_fastq.output.fastq_1),
        fastq_2 = str(rules._igseqr_input_fastq.output.fastq_2),
        fastq_1_real = CFG["inputs"]["sample_fastq_1"],  # Prevent premature deletion of fastqs marked as temp
        fastq_2_real = CFG["inputs"]["sample_fastq_2"],
        hisat_ref_ready = str(rules._igseqr_get_hisat_ref.output.complete),
    output:
        _IGSEQR_RUN_OUTPUTS,
    log:
        stdout = CFG["logs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/igseqr_run.stdout.log",
        stderr = CFG["logs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}/igseqr_run.stderr.log",
    params:
        opts      = CFG["options"]["igseqr_run"],
        out_dir   = CFG["dirs"]["igseqr"] + "{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}",
        hisat_ref = HISAT_REF_DIR + "/" + HISAT_REF_VERSION + "/genome",
        chain_arg = CHAIN_ARG,
    resources:
        **CFG["resources"]["igseqr_run"]
    conda:
        CFG["conda_envs"]["igseqr"]
    container:
        CFG["container_envs"]["igseqr"]
    threads:
        CFG["threads"]["igseqr_run"]
    shell:
        op.as_one_line("""
        igseqr
        --fastq1 {input.fastq_1}
        --fastq2 {input.fastq_2}
        --hisat_ref {params.hisat_ref}
        --out {params.out_dir}
        --sample {wildcards.sample_id}
        --cores {threads}
        --chain {params.chain_arg}
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _igseqr_output_files:
    input:
        _IGSEQR_RUN_OUTPUTS,
    output:
        _IGSEQR_OUTPUT_FILES,
    run:
        for src, dst in zip(input, output):
            op.relative_symlink(src, dst, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _igseqr_all:
    input:
        expand(
            expand(
                [
                    CFG["dirs"]["outputs"] + "fasta/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_transcripts.fasta",
                    CFG["dirs"]["outputs"] + "tsv/{seq_type}--" + HISAT_REF_VERSION + "/{sample_id}_{chain}_report.tsv",
                ],
                chain=CHAINS,
                allow_missing=True,
            ),
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
