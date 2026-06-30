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

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cnaqc"]`
CFG = op.setup_module(
    name = "cnaqc",
    version = "1.0",
    subdirectories = ["inputs", "cnaqc", "outputs"],
)

# Map lcr-modules genome build names to the strings CNAqc's init(ref=) accepts.
# CNAqc accepts: "hg19", "GRCh37", "hg38", "GRCh38"
CNAQC_GENOME_MAP = {
    "hg19":          "hg19",
    "grch37":        "GRCh37",
    "hs37d5":        "GRCh37",
    "hg38":          "hg38",
    "grch38":        "GRCh38",
    "grch38-legacy": "GRCh38",
}

# Define rules to be run locally when using a compute cluster
localrules:
    _cnaqc_input_subclones,
    _cnaqc_input_cellularity_ploidy,
    _cnaqc_input_maf,
    _cnaqc_output,
    _cnaqc_all


##### RULES #####


# Symlink subclones file into module inputs
rule _cnaqc_input_subclones:
    input:
        sub = CFG["inputs"]["subclones"]
    output:
        sub = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_subclones.txt"
    run:
        op.absolute_symlink(input.sub, output.sub)


# Symlink cellularity_ploidy file into module inputs
rule _cnaqc_input_cellularity_ploidy:
    input:
        cp = CFG["inputs"]["cellularity_ploidy"]
    output:
        cp = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cellularity_ploidy.txt"
    run:
        op.absolute_symlink(input.cp, output.cp)


# Symlink MAF file into module inputs
rule _cnaqc_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Run CNAqc analysis
rule _cnaqc_run:
    input:
        sub = str(rules._cnaqc_input_subclones.output.sub),
        cp  = str(rules._cnaqc_input_cellularity_ploidy.output.cp),
        maf = str(rules._cnaqc_input_maf.output.maf),
        script = CFG["inputs"]["cnaqc_script"]
    output:
        plot    = CFG["dirs"]["cnaqc"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc_qc.pdf",
        metrics = CFG["dirs"]["cnaqc"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc_metrics.tsv"
    log:
        stdout = CFG["logs"]["cnaqc"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc.stdout.log",
        stderr = CFG["logs"]["cnaqc"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc.stderr.log"
    params:
        ref                  = lambda w: CNAQC_GENOME_MAP[w.genome_build],
        min_depth            = CFG["options"]["min_depth"],
        min_muts_per_segment = CFG["options"]["min_muts_per_segment"],
        run_peak_analysis    = CFG["options"]["run_peak_analysis"]
    conda:
        CFG["conda_envs"]["cnaqc"]
    container:
        CFG["container_envs"]["cnaqc"]
    threads:
        CFG["threads"]["cnaqc"]
    resources:
        **CFG["resources"]["cnaqc"]
    shell:
        """
        Rscript --vanilla {input.script} \
            --subclones {input.sub} \
            --cellularity_ploidy {input.cp} \
            --maf {input.maf} \
            --tumour_id {wildcards.tumour_id} \
            --out_plot {output.plot} \
            --out_metrics {output.metrics} \
            --ref {params.ref} \
            --min_depth {params.min_depth} \
            --min_muts_per_segment {params.min_muts_per_segment} \
            --run_peak_analysis {params.run_peak_analysis} \
            > {log.stdout} 2> {log.stderr}
        """


# Symlink final output files into 99-outputs/
rule _cnaqc_output:
    input:
        plot    = str(rules._cnaqc_run.output.plot),
        metrics = str(rules._cnaqc_run.output.metrics)
    output:
        plot    = CFG["dirs"]["outputs"] + "pdf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc_qc.pdf",
        metrics = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_cnaqc_metrics.tsv"
    run:
        op.relative_symlink(input.plot, output.plot, in_module=True)
        op.relative_symlink(input.metrics, output.metrics, in_module=True)


# Generates the target sentinels for each run, which generate the symlinks
rule _cnaqc_all:
    input:
        expand(
            [
                str(rules._cnaqc_output.output.plot),
                str(rules._cnaqc_output.output.metrics),
            ],
            zip,
            seq_type   = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id  = CFG["runs"]["tumour_sample_id"],
            normal_id  = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
