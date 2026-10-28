#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


import oncopipe as op

min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
        '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
        "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
    )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

CFG = op.setup_module(
    name = "maf_diff",
    version = "1.0",
    subdirectories = ["inputs", "maf_diff", "outputs"],
)

CALLER1 = CFG["options"]["caller1_name"]
CALLER2 = CFG["options"]["caller2_name"]

localrules:
    _maf_diff_input_maf1,
    _maf_diff_input_maf2,
    _maf_diff_output_maf,
    _maf_diff_output_stats,
    _maf_diff_all


##### RULES #####


rule _maf_diff_input_maf1:
    input:
        maf = CFG["inputs"]["maf1"]
    output:
        maf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf1.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _maf_diff_input_maf2:
    input:
        maf = CFG["inputs"]["maf2"]
    output:
        maf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf2.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _maf_diff_run:
    input:
        maf1 = str(rules._maf_diff_input_maf1.output.maf),
        maf2 = str(rules._maf_diff_input_maf2.output.maf)
    output:
        caller1_only = CFG["dirs"]["maf_diff"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + CALLER1 + "-only.maf",
        caller2_only = CFG["dirs"]["maf_diff"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + CALLER2 + "-only.maf",
        stats = CFG["dirs"]["maf_diff"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stats.tsv"
    log:
        CFG["logs"]["maf_diff"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/maf_diff.log"
    conda:
        CFG["conda_envs"]["pandas"]
    threads:
        CFG["threads"]["maf_diff"]
    resources:
        **CFG["resources"]["maf_diff"]
    params:
        caller1_name = CALLER1,
        caller2_name = CALLER2,
        key_cols = CFG["options"]["key_cols"],
        driver_genes = CFG["options"]["driver_genes"],
        driver_gene_col = CFG["options"]["driver_gene_col"],
        driver_col = CFG["options"]["driver_col"],
        driver_col_value = CFG["options"]["driver_col_value"],
        maf_gene_col = CFG["options"]["maf_gene_col"],
        maf_vc_col = CFG["options"]["maf_vc_col"],
        coding_variant_classes = CFG["options"]["coding_variant_classes"]
    script:
        "src/maf_diff.py"


rule _maf_diff_output_maf:
    input:
        caller1_only = str(rules._maf_diff_run.output.caller1_only),
        caller2_only = str(rules._maf_diff_run.output.caller2_only)
    output:
        caller1_only = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + CALLER1 + "-only.maf",
        caller2_only = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}." + CALLER2 + "-only.maf"
    run:
        op.relative_symlink(input.caller1_only, output.caller1_only, in_module=True)
        op.relative_symlink(input.caller2_only, output.caller2_only, in_module=True)


rule _maf_diff_output_stats:
    input:
        stats = str(rules._maf_diff_run.output.stats)
    output:
        stats = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stats.tsv"
    run:
        op.relative_symlink(input.stats, output.stats, in_module=True)


rule _maf_diff_all:
    input:
        expand(
            [
                str(rules._maf_diff_output_maf.output.caller1_only),
                str(rules._maf_diff_output_maf.output.caller2_only),
                str(rules._maf_diff_output_stats.output.stats),
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        )


##### CLEANUP #####


op.cleanup_module(CFG)
