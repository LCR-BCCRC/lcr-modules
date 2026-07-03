#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Houman Layegh Mirhosseini
# Module Author:    Houman Layegh Mirhosseini
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import os
import re
import glob
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import sys
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
                '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
                "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
                )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section


# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["drivemutr"]`
#
# NOTE: DriveMuTR is a *cohort-level, per-gene* module rather than a per-sample
# or per-tumour/normal module. It consumes a single MAF + copy-number matrix
# spanning the whole cohort, scatters one job per target gene via a checkpoint,
# and gathers the per-gene results into UCSC custom tracks. Its working
# wildcards are therefore `{gene}` (produced by the checkpoint) and `{lam}`
# (one hierarchical-grouping lambda value), not the LCR standard
# `seq_type`/`genome_build`/`sample_id`. `genome_build`/`pathology` are supplied
# through `reference_params` in the config.
CFG = op.setup_module(
    name = "drivemutr",
    version = "1.0",
    subdirectories = ["inputs", "io", "cadd", "per_gene", "ucsc_tracks", "outputs"],
)

# Shorthand for the checkpoint scatter directory
PER_GENE = CFG["dirs"]["per_gene"]

# Tracks emitted by the aggregation step (step 15)
_UCSC_TRACKS = ["Mutation_Points", "Mutation_Blocks", "Transcription_Factors"]


def _lambda_suffix(v):
    """1 -> 'lambda_1', 1.5 -> 'lambda_1_5'."""
    return "lambda_" + str(v).replace(".", "_")

def _lambda_suffixes():
    """The list of lambda suffixes from CFG['options']['lambda']."""
    return [_lambda_suffix(v) for v in CFG["options"]["lambda"]]

def per_gene_data(wildcards):
    """Collect the terminal per-gene RDS files produced after the checkpoint."""
    checkpoint_dir = checkpoints._drivemutr_merge_and_split_genes.get(**wildcards).output[0]
    genes = glob_wildcards(os.path.join(checkpoint_dir, "{gene}/data.rds")).gene
    return expand(os.path.join(checkpoint_dir, "{gene}", "final_results.rds"), gene=genes)

def _prev_step_rds(wildcards):
    """Return the output of step 4 or 5 depending on custom_sliding_window."""
    stem = "sliding_window" if CFG["options"]["custom_sliding_window"] else "lambda_annotated"
    return PER_GENE + wildcards.gene + "/" + stem + ".rds"


wildcard_constraints:
    gene = r"[^/]+",
    track = "|".join([re.escape(t) for t in _UCSC_TRACKS]),
    lam = r"lambda_[0-9A-Za-z_]+",


# Define rules to be run locally when using a compute cluster
localrules:
    _drivemutr_input,
    _drivemutr_output_tracks,
    _drivemutr_all,


##### RULES #####


# Symlinks the cohort-level input files into the module results directory (under '00-inputs/').
# Optional inputs (sample_id_map, cadd_dir) are passed straight from the config to the rules
# that use them, mirroring the battenberg `reference_path` pattern.
rule _drivemutr_input:
    input:
        ssm                  = CFG["inputs"]["ssm_maf"],
        cnv                  = CFG["inputs"]["cnv_matrix"],
        coexpression_modules = CFG["inputs"]["wgcna_coexpression_modules"],
        filtered_expression  = CFG["inputs"]["wgcna_filtered_expression"],
        tf_names             = CFG["inputs"]["tf_names_file"],
    output:
        ssm                  = CFG["dirs"]["inputs"] + "ssm/mutations.maf",
        cnv                  = CFG["dirs"]["inputs"] + "cnv/cn_matrix.tsv",
        coexpression_modules = CFG["dirs"]["inputs"] + "wgcna/coexpression_modules.rds",
        filtered_expression  = CFG["dirs"]["inputs"] + "wgcna/filtered_expression.tsv",
        tf_names             = CFG["dirs"]["inputs"] + "tf/tf_names.txt",
    run:
        op.absolute_symlink(input.ssm, output.ssm)
        op.absolute_symlink(input.cnv, output.cnv)
        op.absolute_symlink(input.coexpression_modules, output.coexpression_modules)
        op.absolute_symlink(input.filtered_expression, output.filtered_expression)
        op.absolute_symlink(input.tf_names, output.tf_names)


# Step 1: retrieve / format input data (collect SSM + CN in aSHM / custom target regions)
rule _drivemutr_annotate_io:
    input:
        ssm = str(rules._drivemutr_input.output.ssm),
        cnv = str(rules._drivemutr_input.output.cnv),
    output:
        rds = temp(CFG["dirs"]["io"] + "target_regions.rds"),
    log:
        CFG["logs"]["io"] + "annotate_io.log",
    params:
        pathology          = CFG["reference_params"]["pathology"],
        genome_build       = CFG["reference_params"]["genome_build"],
        genes_regions_list = CFG["reference_params"]["genes_regions_list"],
        ashm_regions       = CFG["reference_params"]["ashm_regions"],
        sample_id_map_path = CFG["inputs"]["sample_id_map"] or "",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["annotate_io"]
    resources:
        **CFG["resources"]["annotate_io"]
    script:
        "src/run_annotate_drivemutr_io.R"


# Step 2: annotate mutations with CADD scores
rule _drivemutr_get_cadd_scores:
    input:
        rds = str(rules._drivemutr_annotate_io.output.rds),
    output:
        rds = temp(CFG["dirs"]["cadd"] + "annotated_target_regions.rds"),
    log:
        CFG["logs"]["cadd"] + "get_cadd_scores.log",
    params:
        cadd_path = CFG["inputs"]["cadd_dir"] or "",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["get_cadd_scores"]
    resources:
        **CFG["resources"]["get_cadd_scores"]
    script:
        "src/run_get_annotate_cadd_scores.R"


# Step 3: merge multi-region rows per gene and scatter one job per gene (checkpoint)
checkpoint _drivemutr_merge_and_split_genes:
    input:
        rds = str(rules._drivemutr_get_cadd_scores.output.rds),
    output:
        directory(PER_GENE),
    log:
        CFG["logs"]["per_gene"] + "merge_and_split_genes.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["merge_and_split_genes"]
    resources:
        **CFG["resources"]["merge_and_split_genes"]
    script:
        "src/run_merge_multi_ashm_regions_gene.R"


# Step 4: hierarchical grouping of mutation positions (per gene)
rule _drivemutr_annotate_lambda:
    input:
        rds = PER_GENE + "{gene}/data.rds",
    output:
        rds             = temp(PER_GENE + "{gene}/lambda_annotated.rds"),
        height_plot_pdf = PER_GENE + "{gene}/height_plot.pdf",
    log:
        CFG["logs"]["per_gene"] + "{gene}/annotate_lambda.log",
    params:
        lam = CFG["options"]["lambda"],
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["annotate_lambda"]
    resources:
        **CFG["resources"]["annotate_lambda"]
    script:
        "src/run_annotate_lambda_hierarchical_grouping.R"


# Step 5 (optional): sliding-window sub-grouping (per gene)
if CFG["options"]["custom_sliding_window"]:
    rule _drivemutr_annotate_sliding_windows:
        input:
            rds = str(rules._drivemutr_annotate_lambda.output.rds),
        output:
            rds = temp(PER_GENE + "{gene}/sliding_window.rds"),
        log:
            CFG["logs"]["per_gene"] + "{gene}/sliding_window.log",
        params:
            window_size         = CFG["options"]["sliding_window"]["window_size"],
            max_range           = CFG["options"]["sliding_window"]["max_range"],
            slide_by            = CFG["options"]["sliding_window"]["slide_by"],
            min_mutations_count = CFG["options"]["sliding_window"]["min_mutations_count"],
            min_group_range     = CFG["options"]["sliding_window"]["min_group_range"],
        conda:
            CFG["conda_envs"]["drivemutr"]
        container:
            CFG["container_envs"]["drivemutr"]
        threads:
            CFG["threads"]["annotate_sliding_windows"]
        resources:
            **CFG["resources"]["annotate_sliding_windows"]
        script:
            "src/run_annotate_sliding_windows.R"


# Step 6: adjust CADD scores per lambda group (per gene)
rule _drivemutr_adjust_cadd:
    input:
        rds = _prev_step_rds,
    output:
        rds = temp(PER_GENE + "{gene}/cadd_adjusted.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/adjust_cadd.log",
    params:
        lam = CFG["options"]["lambda"],
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["adjust_cadd"]
    resources:
        **CFG["resources"]["adjust_cadd"]
    script:
        "src/run_adjust_cadd.R"


# Step 7: assign WGCNA module colour and expression (per gene)
rule _drivemutr_assign_module:
    input:
        gene_data            = str(rules._drivemutr_adjust_cadd.output.rds),
        coexpression_modules = str(rules._drivemutr_input.output.coexpression_modules),
        filtered_expression  = str(rules._drivemutr_input.output.filtered_expression),
    output:
        rds = temp(PER_GENE + "{gene}/module_annotated.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/assign_module.log",
    params:
        sample_id_map_path = CFG["inputs"]["sample_id_map"] or "",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["assign_module"]
    resources:
        **CFG["resources"]["assign_module"]
    script:
        "src/run_assign_module_color_and_expression.R"


# Step 8: filter to matched DNA + RNA samples (per gene)
rule _drivemutr_filter_matched_dna_rna:
    input:
        rds = str(rules._drivemutr_assign_module.output.rds),
    output:
        rds = temp(PER_GENE + "{gene}/dna_rna_filtered.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/filter_matched_dna_rna.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["filter_matched_dna_rna"]
    resources:
        **CFG["resources"]["filter_matched_dna_rna"]
    script:
        "src/run_filter_matched_DNA_RNA_samples.R"


# Step 9: build mutation foci matrix (per gene)
rule _drivemutr_build_foci_matrix:
    input:
        rds = str(rules._drivemutr_filter_matched_dna_rna.output.rds),
    output:
        rds = temp(PER_GENE + "{gene}/foci_matrix.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/build_mutation_foci_matrix.log",
    params:
        lam = CFG["options"]["lambda"],
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["build_foci_matrix"]
    resources:
        **CFG["resources"]["build_foci_matrix"]
    script:
        "src/run_build_mutation_foci_matrix.R"


# Step 10: module eigengene + linear model residuals (per gene)
rule _drivemutr_module_model:
    input:
        rds = str(rules._drivemutr_build_foci_matrix.output.rds),
    output:
        rds = temp(PER_GENE + "{gene}/module_model.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/module_model.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["module_model"]
    resources:
        **CFG["resources"]["module_model"]
    script:
        "src/run_module_model.R"


# Step 11: genomic model, significant foci, and coefficient annotation (per gene)
rule _drivemutr_genomic_models:
    input:
        rds = str(rules._drivemutr_module_model.output.rds),
    output:
        rds = temp(PER_GENE + "{gene}/genomic_models.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/genomic_models.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["genomic_models"]
    resources:
        **CFG["resources"]["genomic_models"]
    script:
        "src/run_genomic_models.R"


# Step 12: RuleFit model and Shapley interaction analysis (per gene)
rule _drivemutr_rulefit_model:
    input:
        rds = str(rules._drivemutr_genomic_models.output.rds),
    output:
        rds      = temp(PER_GENE + "{gene}/rulefit_model.rds"),
        shap_pdf = PER_GENE + "{gene}/shap_plots.pdf",
    log:
        CFG["logs"]["per_gene"] + "{gene}/rulefit_model.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["rulefit_model"]
    resources:
        **CFG["resources"]["rulefit_model"]
    script:
        "src/run_rulefit_model.R"


# Step 13: TF motif analysis and Fisher enrichment per significant focus (per gene)
rule _drivemutr_get_tf_results:
    input:
        gene_data           = str(rules._drivemutr_rulefit_model.output.rds),
        filtered_expression = str(rules._drivemutr_input.output.filtered_expression),
        tf_names_file       = str(rules._drivemutr_input.output.tf_names),
    output:
        rds = temp(PER_GENE + "{gene}/tf_results.rds"),
    log:
        CFG["logs"]["per_gene"] + "{gene}/get_tf_results.log",
    params:
        expression_threshold = CFG["options"]["tf_expression_threshold"],
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["get_tf_results"]
    resources:
        **CFG["resources"]["get_tf_results"]
    script:
        "src/run_get_tf_results.R"


# Step 14: sanity-check boxplot - gene expression by mutation foci group (per gene)
rule _drivemutr_sanity_check_plot:
    input:
        rds = str(rules._drivemutr_get_tf_results.output.rds),
    output:
        rds              = PER_GENE + "{gene}/final_results.rds",
        sanity_check_pdf = PER_GENE + "{gene}/sanity_check_plot.pdf",
    log:
        CFG["logs"]["per_gene"] + "{gene}/sanity_check_plot.log",
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["sanity_check_plot"]
    resources:
        **CFG["resources"]["sanity_check_plot"]
    script:
        "src/run_annotate_sanity_check_plot.R"


# Step 15: aggregate all per-gene results into UCSC custom track TSV files (gather over genes)
rule _drivemutr_get_ucsc_tracks:
    input:
        per_gene_data,
    output:
        expand(
            CFG["dirs"]["ucsc_tracks"] + "{track}_{lam}.tsv",
            track = _UCSC_TRACKS,
            lam = _lambda_suffixes(),
        ),
    log:
        CFG["logs"]["ucsc_tracks"] + "get_ucsc_tracks.log",
    params:
        outdir = CFG["dirs"]["ucsc_tracks"],
    conda:
        CFG["conda_envs"]["drivemutr"]
    container:
        CFG["container_envs"]["drivemutr"]
    threads:
        CFG["threads"]["get_ucsc_tracks"]
    resources:
        **CFG["resources"]["get_ucsc_tracks"]
    script:
        "src/run_get_ucsc_custom_tracks.R"


# Symlinks the final UCSC custom tracks into the module results directory (under '99-outputs/')
rule _drivemutr_output_tracks:
    input:
        tsv = CFG["dirs"]["ucsc_tracks"] + "{track}_{lam}.tsv",
    output:
        tsv = CFG["dirs"]["outputs"] + "ucsc_custom_tracks/{track}_{lam}.tsv",
    run:
        op.relative_symlink(input.tsv, output.tsv)


# Generates the target sentinels for the module
rule _drivemutr_all:
    input:
        expand(
            str(rules._drivemutr_output_tracks.output.tsv),
            track = _UCSC_TRACKS,
            lam = _lambda_suffixes(),
        ),


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
