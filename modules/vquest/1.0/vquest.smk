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
# `CFG` is a shortcut to `config["lcr-modules"]["vquest"]`
CFG = op.setup_module(
    name = "vquest",
    version = "1.0",
    subdirectories = ["inputs", "vquest", "outputs"],
)

CHAINS = CFG["options"]["chains"]
if isinstance(CHAINS, str):
    CHAINS = CHAINS.split()

VALID_CHAINS = ["IGH", "IGK", "IGL", "IGKL", "TRA", "TRB", "TRD", "TRG"]
assert all(chain in VALID_CHAINS for chain in CHAINS), (
    "Config 'chains' must be a list from: IGH, IGK, IGL, IGKL, TRA, TRB, TRD, TRG. "
    "Use 'IGKL' for combined kappa/lambda FASTA output from igseqr."
)

GLYCO_SOURCE = CFG["options"]["glycosylation_source"]
GLYCO_AA_TSV_PATTERN = CFG["inputs"]["sample_source_tsv"] if GLYCO_SOURCE == "mixcr_tsv" else None

# Maps chain wildcard to the IMGT V-QUEST receptorOrLocusType parameter.
receptor_type_dict = {
    "IGH": "IG", "IGK": "IG", "IGL": "IG", "IGKL": "IG",
    "TRA": "TR", "TRB": "TR", "TRD": "TR",  "TRG": "TR",
}

# Define rules to be run locally when using a compute cluster
localrules:
    _vquest_input_fasta,
    _vquest_output_tsv,
    _vquest_annotated_tsv,
    _vquest_output_annotated_tsv,
    _vquest_merge_tsv,
    _vquest_output_merged_tsv,
    _vquest_merge_annotated_tsv,
    _vquest_output_merged_annotated_tsv,
    _vquest_all,


##### RULES #####


# Symlinks the input FASTA into the module results directory (under '00-inputs/')
rule _vquest_input_fasta:
    input:
        fasta = CFG["inputs"]["sample_fasta"]
    output:
        fasta = CFG["dirs"]["inputs"] + "fasta/{seq_type}/{sample_id}.{chain}.fasta"
    run:
        op.absolute_symlink(input.fasta, output.fasta)


# Submits sequences to IMGT V-QUEST and writes the AIRR TSV output.
# NOTE: This rule requires outbound HTTPS access to www.imgt.org from the
# executing node. On clusters where compute nodes lack internet access,
# run with --cluster-config pointing jobs to a login/submit node, or use
# a workflow proxy. The vquest package batches sequences automatically
# (50 per request) so arbitrarily large FASTA files are handled correctly.
rule _vquest_run:
    input:
        fasta  = str(rules._vquest_input_fasta.output.fasta),
        script = CFG["scripts"]["run_vquest"],
    output:
        tsv = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.vquest_airr.tsv"
    log:
        stdout = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/vquest_run.stdout.log",
        stderr = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/vquest_run.stderr.log",
    params:
        species       = CFG["options"]["species"],
        receptor_type = lambda wildcards: receptor_type_dict[wildcards.chain],
        molecule_type = CFG["options"]["molecule_type"],
    wildcard_constraints:
        chain = "|".join(CHAINS),
    threads:
        CFG["threads"]["vquest_run"]
    resources:
        **CFG["resources"]["vquest_run"]
    conda:
        CFG["conda_envs"]["vquest"]
    shell:
        op.as_one_line("""
        python {input.script}
        --fasta {input.fasta}
        --species {params.species}
        --receptor_type {params.receptor_type}
        --molecule_type {params.molecule_type}
        --output {output.tsv}
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the AIRR TSV into the module results directory (under '99-outputs/')
rule _vquest_output_tsv:
    input:
        tsv = str(rules._vquest_run.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.vquest_airr.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module=True)


# Translates input FASTA and annotates acquired N-linked glycosylation sites
# (NxS/T, x != Pro) with IMGT unique Lefranc numbering via ANARCI.
# In MiXCR mode (glycosylation_source: "mixcr_tsv"), AA sequences are assembled
# from the aaSeq* fields in sample_source_tsv instead of translating the FASTA.
rule _vquest_annotate_glycosylation:
    input:
        fasta  = str(rules._vquest_input_fasta.output.fasta),
        script = CFG["scripts"]["annotate_glycosylation"],
        source_tsv = lambda wildcards: (
            GLYCO_AA_TSV_PATTERN.format(
                seq_type=wildcards.seq_type,
                sample_id=wildcards.sample_id,
                chain=wildcards.chain,
            ) if GLYCO_SOURCE == "mixcr_tsv" else []
        ),
    output:
        tsv = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.glycosylation.tsv"
    log:
        stdout = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/annotate_glycosylation.stdout.log",
        stderr = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/annotate_glycosylation.stderr.log",
    params:
        source_tsv_flag = lambda wildcards, input: (
            f"--source_tsv {input.source_tsv}" if GLYCO_SOURCE == "mixcr_tsv" else ""
        ),
    wildcard_constraints:
        chain = "|".join(CHAINS),
    threads:
        CFG["threads"]["annotate_glycosylation"]
    resources:
        **CFG["resources"]["annotate_glycosylation"]
    conda:
        CFG["conda_envs"]["glycosylation"]
    shell:
        op.as_one_line("""
        python {input.script}
        --fasta {input.fasta}
        {params.source_tsv_flag}
        --output {output.tsv}
        > {log.stdout} 2> {log.stderr}
        """)

# Merges glycosylation annotation into the V-QUEST AIRR TSV (always produced)
rule _vquest_annotated_tsv:
    input:
        vquest = str(rules._vquest_run.output.tsv),
        glyco  = str(rules._vquest_annotate_glycosylation.output.tsv),
    output:
        annotated = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.vquest_airr.annotated.tsv"
    params:
        script = CFG["scripts"]["merge_tsv"],
    wildcard_constraints:
        chain = "|".join(CHAINS),
    shell:
        op.as_one_line("""
        python {params.script}
        --base {input.vquest}
        --annotation {input.glyco}
        --base_key sequence_id
        --annot_key sequence_id
        --output {output.annotated}
        """)

rule _vquest_output_annotated_tsv:
    input:
        annotated = str(rules._vquest_annotated_tsv.output.annotated)
    output:
        annotated = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.vquest_airr.annotated.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.relative_symlink(input.annotated, output.annotated, in_module=True)


# Merges V-QUEST AIRR annotations into the source TSV
rule _vquest_merge_tsv:
    input:
        annotation = str(rules._vquest_run.output.tsv),
        source_tsv = CFG["inputs"]["sample_source_tsv"],
    output:
        merged = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.vquest_airr.merged.tsv"
    params:
        script     = CFG["scripts"]["merge_tsv"],
        source_key = CFG["options"]["source_id_column"],
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

rule _vquest_output_merged_tsv:
    input:
        merged = str(rules._vquest_merge_tsv.output.merged)
    output:
        merged = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.vquest_airr.merged.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.relative_symlink(input.merged, output.merged, in_module=True)

# Merges source TSV with vquest_airr.annotated (V-QUEST + glyco combined)
rule _vquest_merge_annotated_tsv:
    input:
        annotation = str(rules._vquest_annotated_tsv.output.annotated),
        source_tsv = CFG["inputs"]["sample_source_tsv"],
    output:
        merged = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.vquest_airr.merged.annotated.tsv"
    params:
        script     = CFG["scripts"]["merge_tsv"],
        source_key = CFG["options"]["source_id_column"],
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

rule _vquest_output_merged_annotated_tsv:
    input:
        merged = str(rules._vquest_merge_annotated_tsv.output.merged)
    output:
        merged = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.vquest_airr.merged.annotated.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.relative_symlink(input.merged, output.merged, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _vquest_all:
    input:
        expand(
            expand(
                [
                    str(rules._vquest_output_tsv.output.tsv),
                    str(rules._vquest_output_annotated_tsv.output.annotated),
                    str(rules._vquest_output_merged_tsv.output.merged),
                    str(rules._vquest_output_merged_annotated_tsv.output.merged),
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
