#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Sierra Gillis
# Module Author:    Sierra Gillis
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
# `CFG` is a shortcut to `config["lcr-modules"]["mutationtimer"]`
CFG = op.setup_module(
    name = "mutationtimer",
    version = "1.0",
    subdirectories = ["inputs", "convert2bed", "liftover", "resolve_overlaps", "mutationtimer", "plot", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _mutationtimer_input_bb,
    _mutationtimer_input_cellularity,
    _mutationtimer_input_augmented_maf,
    _mutationtimer_convert2bed,
    _mutationtimer_liftover,
    _mutationtimer_resolve_overlaps,
    _mutationtimer_symlink_output_tsvs,
    _mutationtimer_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutationtimer_input_bb:
    input:
        bb = CFG["inputs"]["battenberg_txt"]
    output:
        bb = CFG["dirs"]["inputs"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_subclones.txt"
    group:
        "input_lift_resolve"
    run:
        op.absolute_symlink(input.bb, output.bb)

rule _mutationtimer_input_cellularity:
    input:
        cellularity = CFG["inputs"]["cellularity_txt"]
    output:
        cellularity = CFG["dirs"]["inputs"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_cellularity_ploidy.txt"
    group:
        "input_lift_resolve"
    run:
        op.absolute_symlink(input.cellularity, output.cellularity)

rule _mutationtimer_input_augmented_maf:
    input:
        maf = CFG["inputs"]["augmented_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "augmented_maf/all_slms-3--{projection}.maf"
    group:
        "input_lift_resolve"
    run:
        op.absolute_symlink(input.maf, output.maf)


# Subsets subclones.txt input to only necessary cols and outputs as bed for liftover
rule _mutationtimer_convert2bed:
    input:
        bb = str(rules._mutationtimer_input_bb.output.bb)
    output:
        bed = CFG["dirs"]["convert2bed"] + "battenberg--{genome_build}/{tumour_id}--{normal_id}_subclones.bed"
    log:
        stderr = CFG["logs"]["mutationtimer"] + "{genome_build}/{tumour_id}--{normal_id}/convert2bed.stderr.log"
    group:
        "input_lift_resolve"
    shell:
        op.as_one_line("""
        awk -F"\t" -v OFS="\t" '{{print $1,$2,$3,$8,$9,$10,$11,$12,$13}}' {input.bb} > {output.bed} 2> {log.stderr}
        """)

# Get chain for liftover based on genome build
def get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")

rule _mutationtimer_liftover:
    input:
        bed = str(rules._mutationtimer_convert2bed.output.bed),
        chain = get_chain
    output:
        lifted = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_id}--{normal_id}_lifted_{chain}.bed",
        unampped = CFG["dirs"]["liftover"] + "from--{genome_build}/{tumour_id}--{normal_id}_lifted_{chain}.unmapped.bed"
    log:
        stderr = CFG["logs"]["mutationtimer"] + "{genome_build}/{tumour_id}--{normal_id}/liftover_{chain}.stderr.log"
    params:
        liftover_script = CFG["options"]["liftover_script"],
        minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    wildcard_constraints:
        chain = "hg38ToHg19|hg19ToHg38"
    group:
        "input_lift_resolve"
    shell:
        op.as_one_line("""
        bash {params.liftover_script}  BED {input.bed}
        {output.lifted} {input.chain}
        YES {params.minmatch}
        2> {log.stderr}
        """)

# Ensures the correct subclones file, naive or lifted, is used as input
def _prepare_resolve_inputs(wildcards):
    if "38" in wildcards.projection:
        genome_list = ["hg38", "grch38", "hg38-nci", "hg38-panea"]
    else:
        genome_list = ["hg19", "grch37", "hs37d5", "grch37-noalt", "hg19-clc", "hg19-reddy"]

    CFG = config["lcr-modules"]["mutationtimer"]
    tbl = CFG["runs"]
    tumor_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.normal_sample_id == wildcards.normal_id)]["tumour_genome_build"].tolist()

    # build and projection "match"
    if str(tumor_genome_build[0]) in genome_list:
        bed = str(rules._mutationtimer_convert2bed.output.bed).replace("{genome_build}", tumor_genome_build[0])
    # build and projection differ, and projection is hg38, means build was hg19/grch37/hs37d5
    elif "38" in wildcards.projection:
        bed = str(rules._mutationtimer_liftover.output.lifted).replace("{genome_build}", tumor_genome_build[0]).replace("{chain}", "hg19ToHg38")
    # build and projection differ, and projection is not hg38, means build was hg38/grch38
    else:
        bed = str(rules._mutationtimer_liftover.output.lifted).replace("{genome_build}", tumor_genome_build[0]).replace("{chain}", "hg38ToHg19")

    return bed

rule _mutationtimer_resolve_overlaps:
    input:
        bed = _prepare_resolve_inputs,
    output:
        bed = CFG["dirs"]["resolve_overlaps"] + "{projection}/{tumour_id}--{normal_id}_resolved.bed",
        removed_in_ties = CFG["dirs"]["resolve_overlaps"] + "{projection}/{tumour_id}--{normal_id}_resolved_removed_in_ties.bed"
    log:
        log = CFG["logs"]["mutationtimer"] + "{projection}/{tumour_id}--{normal_id}/resolve_overlaps.stderr.log"
    params:
        script = CFG["options"]["resolve_overlaps_script"]
    conda:
        CFG["conda_envs"]["mutationtimer"]
    group:
        "input_lift_resolve"
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script}
        {input.bed}
        {output.bed}
        {log.log}
        """)


# Ensures the correct subclones file, naive or lifted, is used as input
def _prepare_mt_inputs(wildcards):
    if "38" in wildcards.projection:
        genome_list = ["hg38", "grch38", "hg38-nci", "hg38-panea"]
    else:
        genome_list = ["hg19", "grch37", "hs37d5", "grch37-noalt", "hg19-clc", "hg19-reddy"]

    CFG = config["lcr-modules"]["mutationtimer"]
    tbl = CFG["runs"]
    tumor_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.normal_sample_id == wildcards.normal_id)]["tumour_genome_build"].tolist()

    # build and projection "match"
    if str(tumor_genome_build[0]) in genome_list:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
    # build and projection differ, and projection is hg38, means build was hg19/grch37/hs37d5
    elif "38" in wildcards.projection:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
    # build and projection differ, and projection is not hg38, means build was hg38/grch38
    else:
        cellularity = str(rules._mutationtimer_input_cellularity.output.cellularity).replace("{genome_build}", tumor_genome_build[0])
    return cellularity



rule  _mutationtimer_run:
    input:
        maf = str(rules._mutationtimer_input_augmented_maf.output.maf),
        bb = str(rules._mutationtimer_resolve_overlaps.output.bed),
        cellularity = _prepare_mt_inputs
    output:
        timed_ssm = CFG["dirs"]["mutationtimer"] + "{projection}/{tumour_id}--{normal_id}/{tumour_id}_timed_ssm.{projection}.tsv",
        timed_cna = CFG["dirs"]["mutationtimer"] + "{projection}/{tumour_id}--{normal_id}/{tumour_id}_timed_cna.{projection}.tsv"
    log:
        log = CFG["logs"]["mutationtimer"] + "{projection}/{tumour_id}--{normal_id}/mutationtimer.stderr.log"
    params:
        script = CFG["options"]["mutationtimer_script"],
        n_bootstrap = CFG["options"]["n_bootstrap"]
    conda:
        CFG["conda_envs"]["mutationtimer"]
    threads:
        CFG["threads"]["mutationtimer"]
    resources:
        **CFG["resources"]["mutationtimer"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script}
        {input.bb}
        {input.cellularity}
        {input.maf}
        {params.n_bootstrap}
        {output.timed_ssm}
        {output.timed_cna}
        {wildcards.tumour_id}
        {wildcards.projection}
        {log.log}
        """)

# Creates plots of timed data
rule  _mutationtimer_plot:
    input:
        timed_ssm = str(rules._mutationtimer_run.output.timed_ssm),
        timed_cna = str(rules._mutationtimer_run.output.timed_cna),
    output:
        full_plot = CFG["dirs"]["plot"] + "{projection}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}.{projection}.full.png",
        min_plot = CFG["dirs"]["plot"] + "{projection}/{tumour_id}--{normal_id}/{tumour_id}--{normal_id}.{projection}.min.png"
    log:
        log = CFG["logs"]["plot"] + "{projection}/{tumour_id}--{normal_id}/plot_mutationtimer.stderr.log"
    params:
        script = CFG["options"]["plot_mutationtimer_script"]
    conda:
        CFG["conda_envs"]["plot_mutationtimer"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.script}
        {input.timed_ssm}
        {input.timed_cna}
        {output.full_plot}
        {output.min_plot}
        {wildcards.tumour_id}
        {wildcards.normal_id}
        {wildcards.projection}
        {log.log}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutationtimer_symlink_output_tsvs:
    input:
        timed_ssm = str(rules._mutationtimer_run.output.timed_ssm),
        timed_cna = str(rules._mutationtimer_run.output.timed_cna),
        full_plot = str(rules._mutationtimer_plot.output.full_plot),
        min_plot = str(rules._mutationtimer_plot.output.min_plot)
    output:
        timed_ssm = CFG["dirs"]["outputs"] + "timed_ssm/{projection}/{tumour_id}--{normal_id}_timed_ssm.{projection}.tsv",
        timed_cna = CFG["dirs"]["outputs"] + "timed_cna/{projection}/{tumour_id}--{normal_id}_timed_cna.{projection}.tsv",
        full_plot = CFG["dirs"]["outputs"] + "plots/{projection}/{tumour_id}--{normal_id}.{projection}.full.png",
        min_plot = CFG["dirs"]["outputs"] + "plots/{projection}/{tumour_id}--{normal_id}.{projection}.min.png"
    run:
        op.relative_symlink(input.timed_ssm, output.timed_ssm, in_module= True)
        op.relative_symlink(input.timed_cna, output.timed_cna, in_module= True)
        op.relative_symlink(input.full_plot, output.full_plot, in_module= True)
        op.relative_symlink(input.min_plot, output.min_plot, in_module= True)



# Generates the target sentinels for each run, which generate the symlinks
rule _mutationtimer_all:
    input:
        expand(
            expand(
            [
                str(rules._mutationtimer_symlink_output_tsvs.output.timed_ssm),
                str(rules._mutationtimer_symlink_output_tsvs.output.timed_cna),
                str(rules._mutationtimer_symlink_output_tsvs.output.full_plot),
                str(rules._mutationtimer_symlink_output_tsvs.output.min_plot)
            ],
            zip,
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            allow_missing=True),
            projection=["grch37", "hg38"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
