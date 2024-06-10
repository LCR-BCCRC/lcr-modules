#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

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
# `CFG` is a shortcut to `config["lcr-modules"]["ecotyper"]`
CFG = op.setup_module(
    name = "ecotyper",
    version = "1.0",
    subdirectories = ["inputs", "ecotyper", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _ecotyper_install,
    _ecotyper_input_matrix,
    _ecotyper_input_annotations,
    _ecotyper_output_assignments,
    _ecotyper_all


##### RULES #####

# Download ecotyper from git repo without cloning the repo itself
# Decompress files into the 00-inputs
# Since ecotyper is not tagged with versions, use a specific commit sha
rule _ecotyper_install:
    output:
        ecotyper_script = CFG["dirs"]["inputs"] + "ecotyper/EcoTyper_recovery_bulk.R", # one of the files in the repo
    params:
        url = "https://github.com/digitalcytometry/ecotyper/archive/" + str(CFG["options"]["ecotyper_version"]) + ".tar.gz",
        folder = CFG["dirs"]["inputs"]
    conda:
        CFG["conda_envs"]["wget"]
    group:
        "preprocessing"
    shell:
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}ecotyper/ --strip-components=1
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
# The gene expression matrix
rule _ecotyper_input_matrix:
    input:
        ge_matrix = CFG["inputs"]["ge_matrix"]
    output:
        ge_matrix = CFG["dirs"]["inputs"] + "gene_expression.tsv"
    group:
        "preprocessing"
    run:
        op.absolute_symlink(input.ge_matrix, output.ge_matrix)

# The annotations file
rule _ecotyper_input_annotations:
    input:
        annotations = CFG["inputs"]["annotations"]
    output:
        annotations = CFG["dirs"]["inputs"] + "annotations.tsv"
    group:
        "preprocessing"
    run:
        op.absolute_symlink(input.annotations, output.annotations)


# Preprocess matrix and annotations to prepare run
rule _ecotyper_create_mapping:
    input:
        ge_matrix = str(rules._ecotyper_input_matrix.output.ge_matrix),
        annotations = str(rules._ecotyper_input_annotations.output.annotations)
    output:
        mapping = temp(CFG["dirs"]["ecotyper"] + "mapping.tsv"),
        ge_matrix = temp(CFG["dirs"]["ecotyper"] + "mapped_ge_matrix.tsv"),
        annotations = temp(CFG["dirs"]["ecotyper"] + "mapped_annotations.tsv")
    log:
        stdout = CFG["logs"]["ecotyper"] + "ecotyper_preprocess.stdout.log"
    conda:
        CFG["conda_envs"]["ecotyper"]
    threads:
        CFG["threads"]["processing"]
    resources:
        **CFG["resources"]["processing"]
    group:
        "preprocessing"
    script:
        "src/R/preprocess.R"


# Run ecotyper recovery
rule _ecotyper_run:
    input:
        ecotyper_script = str(rules._ecotyper_install.output.ecotyper_script),
        ge_matrix = str(rules._ecotyper_create_mapping.output.ge_matrix),
        annotations = str(rules._ecotyper_create_mapping.output.annotations)
    output:
        b_cell_assignments = CFG["dirs"]["ecotyper"] + "mapped_ge_matrix/B.cells/state_assignment.txt",
        b_cell_heatmap = CFG["dirs"]["ecotyper"] + "mapped_ge_matrix/B.cells/state_assignment_heatmap.pdf",
        ecotype_assignments = CFG["dirs"]["ecotyper"] + "mapped_ge_matrix/Ecotypes/ecotype_assignment.txt",
        ecotype_heatmap = CFG["dirs"]["ecotyper"] + "mapped_ge_matrix/Ecotypes/heatmap_assigned_samples_viridis.pdf"
    log:
        stdout = CFG["logs"]["ecotyper"] + "ecotyper_run.stdout.log",
        stderr = CFG["logs"]["ecotyper"] + "ecotyper_run.stderr.log"
    params:
        opts = CFG["options"]["annotation_tracks"],
        out_dir = CFG["dirs"]["ecotyper"]
    conda:
        CFG["conda_envs"]["ecotyper"]
    threads:
        CFG["threads"]["ecotyper"]
    resources:
        **CFG["resources"]["ecotyper"]
    group:
        "run"
    shell:
        op.as_one_line("""
        ECOTYPER_INPUT_MATRIX=$(realpath {input.ge_matrix})
            &&
        ECOTYPER_INPUT_ANNOTATIONS=$(realpath {input.annotations})
            &&
        ECOTYPER_STDOUT_LOG=$(realpath {log.stdout})
            &&
        ECOTYPER_STDERR_LOG=$(realpath {log.stderr})
            &&
        ECOTYPER_OUT_DIR=$(readlink -f {params.out_dir})/
            &&
        if [ -e $(echo $ECOTYPER_OUT_DIR)mapped_ge_matrix ]; then rm -R $(echo $ECOTYPER_OUT_DIR)mapped_ge_matrix; fi
            &&
        cd $(dirname {input.ecotyper_script})
            &&
        Rscript --vanilla
        $(basename {input.ecotyper_script})
        -d Lymphoma
        -m $ECOTYPER_INPUT_MATRIX
        -a $ECOTYPER_INPUT_ANNOTATIONS
        {params.opts}
        -o $ECOTYPER_OUT_DIR
        -t {threads}
        > $ECOTYPER_STDOUT_LOG 2> $ECOTYPER_STDERR_LOG
        """)


# Postprocess outputs to return to the same sample ids
rule _ecotyper_postprocess:
    input:
        mapping = str(rules._ecotyper_create_mapping.output.mapping),
        b_cell_assignments = str(rules._ecotyper_run.output.b_cell_assignments)
    output:
        complete = CFG["dirs"]["ecotyper"] + "mapped_ge_matrix/mapping_complete"
    log:
        stdout = CFG["logs"]["ecotyper"] + "ecotyper_postprocess.stdout.log"
    conda:
        CFG["conda_envs"]["ecotyper"]
    threads:
        CFG["threads"]["processing"]
    resources:
        **CFG["resources"]["processing"]
    group:
        "run"
    script:
        "src/R/postprocess.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _ecotyper_output_assignments:
    input:
        mapping_complete = str(rules._ecotyper_postprocess.output.complete),
        b_cell_assignments = str(rules._ecotyper_run.output.b_cell_assignments),
        b_cell_heatmap = str(rules._ecotyper_run.output.b_cell_heatmap),
        ecotype_assignments = str(rules._ecotyper_run.output.ecotype_assignments),
        ecotype_heatmap = str(rules._ecotyper_run.output.ecotype_heatmap)
    output:
        b_cell_assignments = CFG["dirs"]["outputs"] + "bulk_lymphoma_data/B.cells/state_assignment.txt",
        b_cell_heatmap = CFG["dirs"]["outputs"] + "bulk_lymphoma_data/B.cells/state_assignment_heatmap.pdf",
        ecotype_assignments = CFG["dirs"]["outputs"] + "bulk_lymphoma_data/Ecotypes/ecotype_assignment.txt",
        ecotype_heatmap = CFG["dirs"]["outputs"] + "bulk_lymphoma_data/Ecotypes/heatmap_assigned_samples_viridis.pdf"
    run:
        op.relative_symlink(input.b_cell_assignments, output.b_cell_assignments, in_module= True)
        op.relative_symlink(input.b_cell_heatmap, output.b_cell_heatmap, in_module= True)
        op.relative_symlink(input.ecotype_assignments, output.ecotype_assignments, in_module= True)
        op.relative_symlink(input.ecotype_heatmap, output.ecotype_heatmap, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _ecotyper_all:
    input:
        expand(
            [
                str(rules._ecotyper_output_assignments.output.b_cell_assignments),
                str(rules._ecotyper_output_assignments.output.b_cell_heatmap),
                str(rules._ecotyper_output_assignments.output.ecotype_assignments),
                str(rules._ecotyper_output_assignments.output.ecotype_heatmap)
            ]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
