#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     Kostia Dreval


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "oncodriveclustl",
    version = "1.0",
    subdirectories = ["inputs", "oncodriveclustl", "outputs"],
)

assert type(CFG['include_non_coding'])==bool, (
    "Config value for 'include_non_coding' must be set to Boolean value (True/False)"
    "Default is True"
)

assert CFG['options']['clustl'] is not None, (
    "Config value for OncodriveCLUSTL options cannot be None value, if no additional options are required leave as ''"
)

# Define rules to be run locally when using a compute cluster
localrules:
    _oncodriveclustl_input_maf,
    _oncodriveclustl_sample_set,
    _oncodriveclustl_prep_input,
    _oncodriveclustl_txt,
    _oncodriveclustl_all


##### RULES #####

# Symlink master MAF to input directory, easier to extract files this way

rule _oncodriveclustl_input_maf:
    input:
        master_maf = CFG["inputs"]["master_maf"]
    output:
        master_maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/input.maf"
    run:
        op.absolute_symlink(input.master_maf, output.master_maf)

rule _oncodriveclustl_sample_set:
    input:
        sample_set = CFG["inputs"]["sample_sets"]
    output:
        sample_set = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_set, output.sample_set)

rule _oncodriveclustl_prep_input:
    input:
        maf = expand(str(rules._oncodriveclustl_input_maf.output.master_maf), seq_type = CFG["seq_types"], genome_build = CFG["genome_build"]),
        sample_set = str(rules._oncodriveclustl_sample_set.output.sample_set),
        r_script = CFG["prepare_mafs"]
    output:
        txt = temp(CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/{sample_set}.maf")
    params:
        non_coding = str(CFG["include_non_coding"]).upper()
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    shell:
        op.as_one_line("""
        Rscript {input.r_script} 
        {input.maf} {input.sample_set} 
        "$(dirname $(realpath {output.txt}))"
        {wildcards.sample_set} 
        Oncodrive 
        {params.non_coding}
        """)

rule _oncodriveclustl_run:
    input:
        txt = str(rules._oncodriveclustl_prep_input.output.txt)
    output:
        tsv = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/clusters_results.tsv",
        txt = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/elements_results.txt",
        png = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/quantile_quantile_plot.png"
    log:
        stdout = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/oncodriveclustl.stdout.log",
        stderr = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}/oncodriveclustl.stderr.log"
    params:
        local_path = CFG["reference_files_directory"],
        regions = CFG["options"]["regions_file"][CFG["genome_build"].replace("grch37","hg19").replace("grch38","hg38")],
        build = CFG["genome_build"].replace("grch37","hg19").replace("grch38","hg38"),
        opts = CFG["options"]["clustl"]
    threads:
        CFG["threads"]["clustl"]
    resources:
        **CFG["resources"]["clustl"]
    conda: CFG["conda_envs"]["clustl"]
    shell:
        op.as_one_line("""
        export BGDATA_LOCAL={params.local_path} && 
        oncodriveclustl -i {input.txt} -o "$(dirname $(realpath {output.txt}))" 
        -r {params.regions} -g {params.build} 
        --cores {threads} --qqplot {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodriveclustl_txt:
    input:
        tsv = str(rules._oncodriveclustl_run.output.tsv),
        txt = str(rules._oncodriveclustl_run.output.txt),
        png = str(rules._oncodriveclustl_run.output.png)
    output:
        tsv = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}_clusters_results.tsv",
        txt = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}_elements_results.txt",
        png = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}_quantile_quantile_plot.png"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module=True)
        op.relative_symlink(input.txt, output.txt, in_module=True)
        op.relative_symlink(input.png, output.png, in_module=True)

# Generates the target sentinels for each run, which generate the symlinks
rule _oncodriveclustl_all:
    input:
        expand(
            [
                str(rules._oncodriveclustl_txt.output.txt),
                str(rules._oncodriveclustl_txt.output.tsv),
                str(rules._oncodriveclustl_txt.output.png)
            ],
            zip,
            genome_build=CFG["genome_build"],
            sample_set=CFG["sample_set"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
