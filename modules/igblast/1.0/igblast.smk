#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Manuela Cruz
# Module Author:    Laura Hilton
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
# `CFG` is a shortcut to `config["lcr-modules"]["igblast"]`
CFG = op.setup_module(
    name = "igblast",
    version = "1.0",
    subdirectories = ["inputs", "igblast", "outputs"],
)

CHAINS = CFG["options"]["chains"]
if isinstance(CHAINS, str):
    CHAINS = CHAINS.split(" ")

VALID_CHAINS = ["IGH", "IGK", "IGL", "IGKL", "TRA", "TRB", "TRD", "TRG"]
assert all(chain in VALID_CHAINS for chain in CHAINS), (
    "Config 'chains' must be a list from: IGH, IGK, IGL, IGKL, TRA, TRB, TRD, TRG. "
    "Use 'IGKL' for combined kappa/lambda FASTA output from igseqr."
)

RUN_MERGE = CFG["options"]["run_merge"]
SOURCE_TSV_PATTERN = CFG["inputs"]["sample_source_tsv"] if RUN_MERGE else None
SOURCE_ID_COL = CFG["options"]["source_id_column"]

# Maps chain name to igblastn seqtype and germline DB prefix
receptor_dict = {
    "IGH": "ig", "IGK": "ig", "IGL": "ig", "IGKL": "ig",
    "TRA": "tr", "TRB": "tr", "TRD": "tr", "TRG": "tr",
}
run_dict = {
    "IGH": "Ig", "IGK": "Ig", "IGL": "Ig", "IGKL": "Ig",
    "TRA": "TCR", "TRB": "TCR", "TRD": "TCR", "TRG": "TCR",
}

# Define rules to be run locally when using a compute cluster
localrules:
    _igblast_input_fasta,
    _igblast_parse_tsv,
    _igblast_output_tsv,
    _igblast_all,

if RUN_MERGE:
    localrules:
        _igblast_merge_tsv,
        _igblast_output_merged_tsv,


##### RULES #####


# Symlinks the input FASTA into the module results directory (under '00-inputs/')
# Set sample_fasta in your project config to point to the per-chain FASTA files:
#   From mixcr 1.3: "{MIXCR_OUT}/99-outputs/fasta/{seq_type}/mixcr.{sample_id}.clonotypes.{chain}.VDJseq.fasta"
#   From igseqr 1.0: "{IGSEQR_OUT}/99-outputs/fasta/{seq_type}--{genome_build}/{sample_id}_{chain}_transcripts.fasta"
rule _igblast_input_fasta:
    input:
        fasta = CFG["inputs"]["sample_fasta"]
    output:
        fasta = CFG["dirs"]["inputs"] + "fasta/{seq_type}/{sample_id}.{chain}.fasta"
    run:
        op.absolute_symlink(input.fasta, output.fasta)

rule _igblastn_run:
    input:
        fasta = str(rules._igblast_input_fasta.output.fasta),
        ig_db = reference_files("genomes/no_build/igblast/database/imgt_database.success")
    output:
        fmt7 = CFG["dirs"]["igblast"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.igblastn.fmt7"
    params:
        receptor_type = lambda wildcards: run_dict[wildcards.chain],
        aux = reference_files("downloads/igblast/optional_file/human_gl.aux"),
        gdv = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_v.ndb")).replace(".ndb",""),
        gdj = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_j.ndb")).replace(".ndb",""),
        gdd = lambda wildcards: (reference_files("genomes/no_build/igblast/database/imgt_human_" + receptor_dict[wildcards.chain] + "_d.ndb")).replace(".ndb",""),
        run_flags = CFG["options"]["igblast_run"]["run_flags"],
        form = "7 std btop " + CFG["options"]["igblast_run"]["form"]
    conda:
        CFG["conda_envs"]["igblast"]
    container:
        CFG["container_envs"]["igblast"]
    threads:
        CFG["threads"]["igblastn_run"]
    resources:
        **CFG["resources"]["igblastn_run"]
    shell:
        op.as_one_line("""
        igblastn -query {input.fasta} -out {output.fmt7}
        -ig_seqtype {params.receptor_type} -organism human
        -auxiliary_data {params.aux}
        -germline_db_V {params.gdv} -germline_db_J {params.gdj} -germline_db_D {params.gdd}
        {params.run_flags} -outfmt '{params.form}' -domain_system imgt
        """)

# Parse igblastn fmt7 into a human-readable TSV
rule _igblast_parse_tsv:
    input:
        fmt7 = str(rules._igblastn_run.output.fmt7),
        script = CFG["scripts"]["igblast2tsv"]
    output:
        tsv = CFG["dirs"]["igblast"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.igblastn.tsv"
    shell:
        "{input.script} -i {input.fmt7} -o {output.tsv}"

# Symlinks the TSV output into the module results directory (under '99-outputs/')
rule _igblast_output_tsv:
    input:
        tsv = str(rules._igblast_parse_tsv.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.igblastn.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module=True)

# Merges igblastn annotations into the source TSV (optional, controlled by run_merge config)
if RUN_MERGE:
    rule _igblast_merge_tsv:
        input:
            annotation = str(rules._igblast_parse_tsv.output.tsv),
            source_tsv = lambda wildcards: SOURCE_TSV_PATTERN.format(
                seq_type=wildcards.seq_type,
                sample_id=wildcards.sample_id,
                chain=wildcards.chain,
            ),
        output:
            merged = CFG["dirs"]["igblast"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.igblastn.merged.tsv"
        params:
            script     = CFG["scripts"]["merge_tsv"],
            source_key = SOURCE_ID_COL,
        wildcard_constraints:
            chain = "|".join(CHAINS),
        shell:
            op.as_one_line("""
            python {params.script}
            --base {input.source_tsv}
            --annotation {input.annotation}
            --base_key {params.source_key}
            --annot_key sequence_id
            --output {output.merged}
            """)

    rule _igblast_output_merged_tsv:
        input:
            merged = str(rules._igblast_merge_tsv.output.merged)
        output:
            merged = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.igblastn.merged.tsv"
        wildcard_constraints:
            chain = "|".join(CHAINS)
        run:
            op.relative_symlink(input.merged, output.merged, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
if RUN_MERGE:
    rule _igblast_all:
        input:
            expand(
                expand(
                    [
                        str(rules._igblast_output_tsv.output.tsv),
                        str(rules._igblast_output_merged_tsv.output.merged),
                    ],
                    zip,
                    seq_type=CFG["samples"]["seq_type"],
                    sample_id=CFG["samples"]["sample_id"],
                    allow_missing=True),
                chain=CHAINS)
else:
    rule _igblast_all:
        input:
            expand(
                expand(
                    [
                        str(rules._igblast_output_tsv.output.tsv),
                    ],
                    zip,
                    seq_type=CFG["samples"]["seq_type"],
                    sample_id=CFG["samples"]["sample_id"],
                    allow_missing=True),
                chain=CHAINS)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
