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


# Maps chain wildcard to the IMGT V-QUEST receptorOrLocusType parameter.
receptor_type_dict = {
    "IGH": "IG", "IGK": "IG", "IGL": "IG", "IGKL": "IG",
    "TRA": "TR", "TRB": "TR", "TRD": "TR",  "TRG": "TR",
}

# Define rules to be run locally when using a compute cluster
localrules:
    _vquest_input_fasta,
    _vquest_input_source_tsv,
    _vquest_filter_long_seqs,
    _vquest_output_tsv,
    _vquest_output_merged_final,
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


# Symlinks the sample source TSV into the module results directory (under '00-inputs/')
rule _vquest_input_source_tsv:
    input:
        tsv = CFG["inputs"]["sample_source_tsv"]
    output:
        tsv = CFG["dirs"]["inputs"] + "tsv/{seq_type}/{sample_id}.{chain}.source.tsv"
    wildcard_constraints:
        chain = "|".join(CHAINS)
    run:
        op.absolute_symlink(input.tsv, output.tsv)


# Drops sequences longer than 10000 bp, which V-QUEST cannot process.
rule _vquest_filter_long_seqs:
    input:
        fasta = str(rules._vquest_input_fasta.output.fasta),
    output:
        fasta = temp(CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.filtered.fasta")
    wildcard_constraints:
        chain = "|".join(CHAINS),
    shell:
        "awk 'BEGIN{{RS=\">\"; FS=\"\\n\"}} NR>1 {{header=$1; seq=\"\"; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)<=10000) printf \">%s\\n%s\\n\", header, seq}}' {input.fasta} > {output.fasta}"


# Submits sequences to IMGT V-QUEST and writes the AIRR TSV output.
# NOTE: This rule requires outbound HTTPS access to www.imgt.org from the
# executing node. On clusters where compute nodes lack internet access,
# run with --cluster-config pointing jobs to a login/submit node, or use
# a workflow proxy. The vquest package batches sequences automatically
# (50 per request) so arbitrarily large FASTA files are handled correctly.
rule _vquest_run:
    input:
        fasta  = str(rules._vquest_filter_long_seqs.output.fasta),
    output:
        tsv = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.vquest_airr.tsv"
    log:
        stdout = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/vquest_run.stdout.log",
        stderr = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/vquest_run.stderr.log",
    params:
        script        = CFG["scripts"]["run_vquest"],
        species       = CFG["options"]["species"],
        receptor_type = lambda wildcards: receptor_type_dict[wildcards.chain],
        molecule_type = CFG["options"]["molecule_type"],
    wildcard_constraints:
        chain = "|".join(CHAINS),
    retries:
        CFG["options"]["vquest_retries"]
    threads:
        CFG["threads"]["vquest_run"]
    resources:
        **CFG["resources"]["vquest_run"]
    conda:
        CFG["conda_envs"]["vquest"]
    shell:
        op.as_one_line("""
        if [ ! -s {input.fasta} ]; then
            touch {output.tsv} ;
        else
            python {params.script}
            --fasta {input.fasta}
            --species {params.species}
            --receptor_type {params.receptor_type}
            --molecule_type {params.molecule_type}
            --output {output.tsv}
            > {log.stdout} 2> {log.stderr} ;
        fi
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


# Annotates acquired N-linked glycosylation sites (NxS/T, x != Pro) with IMGT
# unique Lefranc numbering via ANARCI. AA sequences are read from the
# sequence_alignment_aa field of the V-QUEST AIRR TSV.
rule _vquest_annotate_glycosylation:
    input:
        fasta      = str(rules._vquest_filter_long_seqs.output.fasta),
        source_tsv = str(rules._vquest_run.output.tsv),
    output:
        tsv = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.glycosylation.tsv"
    log:
        stdout = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/annotate_glycosylation.stdout.log",
        stderr = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/annotate_glycosylation.stderr.log",
    params:
        script = CFG["scripts"]["annotate_glycosylation"],
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
        python {params.script}
        --fasta {input.fasta}
        --source_tsv {input.source_tsv}
        --output {output.tsv}
        > {log.stdout} 2> {log.stderr}
        """)

# Merges vquest_airr + glycosylation columns into source TSV in two passes:
#   1. join glyco into vquest_airr (both keyed on sequence_id)
#   2. join the combined annotation into source_tsv (keyed on source_id_column)
rule _vquest_merge_final:
    input:
        source_tsv = str(rules._vquest_input_source_tsv.output.tsv),
        vquest     = str(rules._vquest_run.output.tsv),
        glyco      = str(rules._vquest_annotate_glycosylation.output.tsv),
    output:
        merged = CFG["dirs"]["vquest"] + "{seq_type}/{sample_id}/{sample_id}.{chain}.merged.tsv"
    log:
        stdout = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/merge_final.stdout.log",
        stderr = CFG["logs"]["vquest"] + "{seq_type}/{sample_id}/{chain}/merge_final.stderr.log",
    params:
        script     = CFG["scripts"]["merge_tsv"],
        source_key = CFG["options"]["source_id_column"],
    wildcard_constraints:
        chain = "|".join(CHAINS),
    shell:
        op.as_one_line("""
        if [ ! -s {input.vquest} ] || [ ! -s {input.source_tsv} ]; then
            touch {output.merged} ;
        else
            tmp=$(mktemp --suffix=.tsv) &&
            python {params.script}
            --base {input.vquest}
            --annotation {input.glyco}
            --base_key sequence_id
            --annot_key sequence_id
            --output $tmp
            > {log.stdout} 2> {log.stderr} &&
            python {params.script}
            --base $tmp
            --annotation {input.source_tsv}
            --base_key sequence_id
            --annot_key {params.source_key}
            --sample_id {wildcards.sample_id}
            --output {output.merged}
            >> {log.stdout} 2>> {log.stderr} &&
            rm $tmp ;
        fi
        """)

rule _vquest_output_merged_final:
    input:
        merged = str(rules._vquest_merge_final.output.merged)
    output:
        merged = CFG["dirs"]["outputs"] + "tsv/{seq_type}/{sample_id}.{chain}.merged.tsv"
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
                    str(rules._vquest_output_merged_final.output.merged),
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
