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
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "ecotyper", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _ecotyper_input_tsv,
    _ecotyper_step_2,
    _ecotyper_output_tsv,
    _ecotyper_all,


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
    group:
        "preprocessing"
    shell:
        op.as_one_line("""
        mkdir {params.folder}/ecotyper
            &&
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}/ecotyper/ --strip-components=1
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
        b_cell_assignments = CFG["dirs"]["ecotyper"] + "RecoveryOutput/bulk_lymphoma_data/B.cells/state_assignment.txt",
        b_cell_heatmap = CFG["dirs"]["ecotyper"] + "RecoveryOutput/bulk_lymphoma_data/B.cells/state_assignment_heatmap.pdf",
        ecotype_assignments = CFG["dirs"]["ecotyper"] + "RecoveryOutput/bulk_lymphoma_data/Ecotypes/ecotype_assignment.txt",
        ecotype_heatmap = CFG["dirs"]["ecotyper"] + "RecoveryOutput/bulk_lymphoma_data/Ecotypes/heatmap_assigned_samples_viridis.pdf"
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
    shell:
        op.as_one_line("""
        cd $(realpath $(dirname {input.ecotyper_script}))
            &&
        Rscript
        {input.ecotyper_script}
        -d Lymphoma
        -m {input.ge_matrix}
        -a {input.annotations}
        {params.opts}
        -o $(realpath $(dirname {params.out_dir}))"/RecoveryOutput"
        -t {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Postprocess outputs to return to the same sample ids
rule _ecotyper_postprocess:
    input:
        mapping = str(rules._ecotyper_create_mapping.output.mapping),
        b_cell_assignments = str(rules._ecotyper_run.output.b_cell_assignments)
    output:
        complete = CFG["dirs"]["ecotyper"] + "RecoveryOutput/bulk_lymphoma_data/mapping_complete"
    conda:
        CFG["conda_envs"]["ecotyper"]
    threads:
        CFG["threads"]["processing"]
    resources:
        **CFG["resources"]["processing"]
    script:
        "src/R/postprocess.R"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _ecotyper_output_tsv:
    input:
        tsv = str(rules._ecotyper_step_2.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _ecotyper_all:
    input:
        expand(
            [
                str(rules._ecotyper_output_tsv.output.tsv),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
