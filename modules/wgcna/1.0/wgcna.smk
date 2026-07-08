#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Houman Layegh Mirhosseini
# Module Author:    Houman Layegh Mirhosseini
# Contributors:     N/A


##### SETUP #####


import oncopipe as op

# Check that the oncopipe dependency is up-to-date
min_oncopipe_version = "1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit(
        "The packaging module dependency is missing. Please install it "
        "('pip install packaging') and ensure you are using the most "
        "up-to-date oncopipe version"
    )

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
        '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
        "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. "
        "Please update oncopipe in your environment" + '\x1b[0m'
    )
    sys.exit(
        "Instructions for updating to the current version of oncopipe are available at "
        "https://lcr-modules.readthedocs.io/en/latest/ (use option 2)"
    )

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["wgcna"]`
CFG = op.setup_module(
    name = "wgcna",
    version = "1.0",
    subdirectories = ["inputs", "wgcna", "outputs"],
)

# Preserve CFG reference for use inside lambda functions
_wgcna_CFG = CFG

# Resolve input mode and pathologies
_wgcna_ge_cfg = CFG["gene_expression"]
_wgcna_mode   = _wgcna_ge_cfg["input_mode"]

if _wgcna_mode not in ("raw", "normalized"):
    raise ValueError(
        f"gene_expression.input_mode must be 'raw' or 'normalized', "
        f"got '{_wgcna_mode}'"
    )

if _wgcna_mode == "raw":
    _raw_pathologies  = _wgcna_ge_cfg.get("pathologies")
    WGCNA_PATHOLOGIES = _raw_pathologies if _raw_pathologies else ["all"]
else:
    WGCNA_PATHOLOGIES = ["user_provided"]

wildcard_constraints:
    pathology = "[^/]+"


localrules:
    _wgcna_input_matrix,
    _wgcna_output_modules,
    _wgcna_all,


##### RULES #####


# Symlink the input expression matrix into the module inputs directory
rule _wgcna_input_matrix:
    input:
        matrix = CFG["inputs"]["expression_matrix"],
    output:
        matrix = CFG["dirs"]["inputs"] + "expression_matrix",
    run:
        op.absolute_symlink(input.matrix, output.matrix)


# Rule 1 (raw mode only): DESeq2 VST normalization + limma batch correction
if _wgcna_mode == "raw":
    rule _wgcna_normalize_expression:
        input:
            raw_matrix = str(rules._wgcna_input_matrix.output.matrix),
        output:
            tsv = CFG["dirs"]["wgcna"] + "{pathology}/normalized_expression.tsv",
        log:
            log = CFG["logs"]["wgcna"] + "{pathology}/normalize_expression.log",
        threads: CFG["threads"]["normalize_expression"]
        resources:
            mem_mb = CFG["mem_mb"]["normalize_expression"],
        container:
            CFG["container_envs"]["wgcna"]
        conda:
            CFG["conda_envs"]["wgcna"]
        params:
            opts                  = CFG["options"]["normalize_expression"],
            pathology             = lambda wc: None if wc.pathology == "all" else wc.pathology,
            raw_matrix_path       = str(rules._wgcna_input_matrix.output.matrix),
            failed_qc_path        = _wgcna_ge_cfg.get("failed_qc_path", ""),
            samples_metadata_path = _wgcna_ge_cfg["samples_metadata_path"],
            plots_dir             = lambda wc: _wgcna_CFG["dirs"]["wgcna"] + wc.pathology + "/plots",
            min_reads             = _wgcna_ge_cfg.get("min_reads", 3000000),
            batches               = _wgcna_ge_cfg.get("batches", ["protocol", "ffpe_or_frozen"]),
            cohort_var            = _wgcna_ge_cfg.get("cohort_var", "cohort"),
            bio_var               = _wgcna_ge_cfg.get("bio_var", "pathology"),
            script_dir            = workflow.basedir + "/src",
        script:
            "src/run_normalize_gene_expression.R"


# Helper: resolve which normalized matrix feeds the filter step
def _wgcna_normalized_input(wildcards):
    if _wgcna_mode == "raw":
        return _wgcna_CFG["dirs"]["wgcna"] + wildcards.pathology + "/normalized_expression.tsv"
    return _wgcna_CFG["dirs"]["inputs"] + "expression_matrix"


# Rule 2: keep genes above median and MAD thresholds
rule _wgcna_filter_variance_genes:
    input:
        expression = _wgcna_normalized_input,
    output:
        tsv = CFG["dirs"]["wgcna"] + "{pathology}/filtered_expression.tsv",
    log:
        log = CFG["logs"]["wgcna"] + "{pathology}/filter_variance_genes.log",
    threads: CFG["threads"]["filter_variance_genes"]
    resources:
        mem_mb = CFG["mem_mb"]["filter_variance_genes"],
    container:
        CFG["container_envs"]["wgcna"]
    conda:
        CFG["conda_envs"]["wgcna"]
    params:
        opts             = CFG["options"]["filter_variance_genes"],
        median_threshold = _wgcna_ge_cfg["median_threshold"],
        mad_threshold    = _wgcna_ge_cfg["mad_threshold"],
        plots_dir        = lambda wc: _wgcna_CFG["dirs"]["wgcna"] + wc.pathology + "/plots",
        script_dir       = workflow.basedir + "/src",
    script:
        "src/run_filter_high_variance_genes.R"


# Rule 3: soft-threshold selection and blockwiseModules WGCNA run
rule _wgcna_get_coexpression_modules:
    input:
        expression = str(rules._wgcna_filter_variance_genes.output.tsv),
    output:
        modules_rds = CFG["dirs"]["wgcna"] + "{pathology}/coexpression_modules.rds",
        modules_tsv = CFG["dirs"]["wgcna"] + "{pathology}/coexpression_modules.tsv",
        network_rds = CFG["dirs"]["wgcna"] + "{pathology}/network.rds",
    log:
        log = CFG["logs"]["wgcna"] + "{pathology}/get_coexpression_modules.log",
    benchmark:
        CFG["dirs"]["wgcna"] + "{pathology}/benchmark_get_coexpression_modules.tsv"
    threads: CFG["threads"]["get_coexpression_modules"]
    resources:
        mem_mb = CFG["mem_mb"]["get_coexpression_modules"],
    container:
        CFG["container_envs"]["wgcna"]
    conda:
        CFG["conda_envs"]["wgcna"]
    params:
        opts         = CFG["options"]["get_coexpression_modules"],
        powers       = _wgcna_ge_cfg["wgcna_powers"],
        cor_method   = _wgcna_ge_cfg["cor_method"],
        network_type = _wgcna_ge_cfg["network_type"],
        plots_dir    = lambda wc: _wgcna_CFG["dirs"]["wgcna"] + wc.pathology + "/plots",
        script_dir   = workflow.basedir + "/src",
    script:
        "src/run_get_co_expressed_genes_modules.R"


# Symlink final outputs into the module outputs directory
rule _wgcna_output_modules:
    input:
        modules_rds = str(rules._wgcna_get_coexpression_modules.output.modules_rds),
        modules_tsv = str(rules._wgcna_get_coexpression_modules.output.modules_tsv),
        network_rds = str(rules._wgcna_get_coexpression_modules.output.network_rds),
    output:
        modules_rds = CFG["dirs"]["outputs"] + "{pathology}/coexpression_modules.rds",
        modules_tsv = CFG["dirs"]["outputs"] + "{pathology}/coexpression_modules.tsv",
        network_rds = CFG["dirs"]["outputs"] + "{pathology}/network.rds",
    run:
        op.relative_symlink(input.modules_rds, output.modules_rds, in_module=True)
        op.relative_symlink(input.modules_tsv, output.modules_tsv, in_module=True)
        op.relative_symlink(input.network_rds, output.network_rds, in_module=True)


# Target rule aggregating all final outputs
rule _wgcna_all:
    input:
        expand(
            [
                str(rules._wgcna_output_modules.output.modules_rds),
                str(rules._wgcna_output_modules.output.modules_tsv),
                str(rules._wgcna_output_modules.output.network_rds),
            ],
            pathology=WGCNA_PATHOLOGIES,
        )


##### CLEANUP #####


op.cleanup_module(CFG)
