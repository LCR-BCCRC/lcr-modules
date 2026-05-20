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

# Define rules to be run locally when using a compute cluster
localrules:
    _igseqr_input_bam,
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


# Symlinks the input BAM and index into the module results directory (under '00-inputs/')
rule _igseqr_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"],
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)


# Run IgSeqR to assemble and quantify immunoglobulin transcripts
rule _igseqr_run:
    input:
        bam = rules._igseqr_input_bam.output.bam,
        bam_real = CFG["inputs"]["sample_bam"],  # Prevent premature deletion of temp bam
        hisat_ref_ready = str(rules._igseqr_get_hisat_ref.output.complete),
    output:
        igh_fasta      = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGH_transcripts.fasta",
        igkl_fasta     = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGKL_transcripts.fasta",
        igh_report     = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGH_report.tsv",
        igkl_report    = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGKL_report.tsv",
        igh_dominant   = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGH_dominant_report.tsv",
        igkl_dominant  = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGKL_dominant_report.tsv",
        igh_tpm_fasta  = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGH_TPM_filtered.fasta",
        igkl_tpm_fasta = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/{sample_id}_IGKL_TPM_filtered.fasta",
    log:
        stdout = CFG["logs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/igseqr_run.stdout.log",
        stderr = CFG["logs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}/igseqr_run.stderr.log",
    params:
        opts      = CFG["options"]["igseqr_run"],
        out_dir   = CFG["dirs"]["igseqr"] + "{seq_type}--{genome_build}/{sample_id}/{hisat_ref_version}",
        hisat_ref = HISAT_REF_DIR + "/" + HISAT_REF_VERSION + "/genome",
    wildcard_constraints:
        hisat_ref_version = HISAT_REF_VERSION
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
        --bam {input.bam}
        --hisat_ref {params.hisat_ref}
        --out {params.out_dir}
        --sample {wildcards.sample_id}
        --cores {threads}
        --chain B
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _igseqr_output_files:
    input:
        igh_fasta      = rules._igseqr_run.output.igh_fasta,
        igkl_fasta     = rules._igseqr_run.output.igkl_fasta,
        igh_report     = rules._igseqr_run.output.igh_report,
        igkl_report    = rules._igseqr_run.output.igkl_report,
        igh_dominant   = rules._igseqr_run.output.igh_dominant,
        igkl_dominant  = rules._igseqr_run.output.igkl_dominant,
        igh_tpm_fasta  = rules._igseqr_run.output.igh_tpm_fasta,
        igkl_tpm_fasta = rules._igseqr_run.output.igkl_tpm_fasta,
    output:
        igh_fasta      = CFG["dirs"]["outputs"] + "fasta/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGH_transcripts.fasta",
        igkl_fasta     = CFG["dirs"]["outputs"] + "fasta/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGKL_transcripts.fasta",
        igh_report     = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGH_report.tsv",
        igkl_report    = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGKL_report.tsv",
        igh_dominant   = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGH_dominant_report.tsv",
        igkl_dominant  = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGKL_dominant_report.tsv",
        igh_tpm_fasta  = CFG["dirs"]["outputs"] + "fasta/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGH_TPM_filtered.fasta",
        igkl_tpm_fasta = CFG["dirs"]["outputs"] + "fasta/{seq_type}--{genome_build}/{hisat_ref_version}/{sample_id}_IGKL_TPM_filtered.fasta",
    wildcard_constraints:
        hisat_ref_version = HISAT_REF_VERSION
    run:
        op.relative_symlink(input.igh_fasta,      output.igh_fasta,      in_module=True)
        op.relative_symlink(input.igkl_fasta,     output.igkl_fasta,     in_module=True)
        op.relative_symlink(input.igh_report,     output.igh_report,     in_module=True)
        op.relative_symlink(input.igkl_report,    output.igkl_report,    in_module=True)
        op.relative_symlink(input.igh_dominant,   output.igh_dominant,   in_module=True)
        op.relative_symlink(input.igkl_dominant,  output.igkl_dominant,  in_module=True)
        op.relative_symlink(input.igh_tpm_fasta,  output.igh_tpm_fasta,  in_module=True)
        op.relative_symlink(input.igkl_tpm_fasta, output.igkl_tpm_fasta, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _igseqr_all:
    input:
        expand(
            expand(
                [
                    rules._igseqr_output_files.output.igh_fasta,
                    rules._igseqr_output_files.output.igkl_fasta,
                    rules._igseqr_output_files.output.igh_report,
                    rules._igseqr_output_files.output.igkl_report,
                ],
                zip,  # Run expand() with zip(), not product()
                seq_type=CFG["samples"]["seq_type"],
                genome_build=CFG["samples"]["genome_build"],
                sample_id=CFG["samples"]["sample_id"],
                allow_missing=True),
            hisat_ref_version=HISAT_REF_VERSION)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
