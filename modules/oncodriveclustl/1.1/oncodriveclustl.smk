#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     Kostia Dreval


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
from datetime import datetime
import numpy as np

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "oncodriveclustl",
    version = "1.1",
    subdirectories = ["inputs", "prepare_mafs", "oncodriveclustl", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _oncodriveclustl_input_maf,
    _oncodriveclustl_sample_set,
    _oncodriveclustl_prep_input,
    _oncodriveclustl_blacklist,
    _oncodriveclustl_format_input,
    _oncodriveclustl_out,
    _oncodriveclustl_genomic_coordinates_out,
    _oncodriveclustl_aggregate,
    _oncodriveclustl_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG["launch_date"]
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later
PREPARE_MAFS = os.path.abspath(config["lcr-modules"]["oncodriveclustl"]["maf_processing"]["prepare_mafs"])

# Symlink master MAF to input directory, easier to extract files this way

rule _oncodriveclustl_input_maf:
    input:
        maf = CFG["inputs"]["input_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "master_maf/{seq_type}--{genome_build}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

rule _oncodriveclustl_sample_set:
    input:
        subsetting_categories = ancient(CFG["inputs"]["subsetting_categories"])
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

checkpoint _oncodriveclustl_prep_input:
    input:
        maf = expand(
            str(rules._oncodriveclustl_input_maf.output.maf),
            seq_type = CFG["samples"]["seq_type"].unique(),
            allow_missing=True),
        subsetting_categories = ancient(str(rules._oncodriveclustl_sample_set.output.subsetting_categories))
    output:
        CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/done"
    log:
        stdout = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/prep_input/prep_input_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["maf_processing"]["include_non_coding"]).upper(),
        mode = "OncodriveCLUSTL",
        metadata_cols = CFG["samples"],
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

rule _oncodriveclustl_blacklist:
    input:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.maf",
        blacklists = CFG["maf_processing"]["blacklists"],
        deblacklist_script = CFG["scripts"]["deblacklist_script"]
    output:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.deblacklisted.maf"
    params:
        drop_threshold = CFG["maf_processing"]["blacklist_drop_threshold"]
    log:
        stdout = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stdout.log",
        stderr = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stderr.log"
    shell:
        op.as_one_line("""
        {input.deblacklist_script}
        --input {input.maf}
        --output {output.maf}
        --drop-threshold {params.drop_threshold}
        --blacklists {input.blacklists}
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodriveclustl_format_input:
    input:
        maf = str(rules._oncodriveclustl_blacklist.output.maf)
    output:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.clustl_input.maf"
    params:
        columns = CFG["format_clustl_input"]["input_columns"],
        additional_commands = CFG["format_clustl_input"]["additional_commands"]
    shell:
        op.as_one_line("""
        cut -f {params.columns} {input.maf} |
        sed 's/Chromosome/CHROMOSOME/' |
        sed 's/Start_Position/POSITION/' |
        sed 's/Reference_Allele/REF/' |
        sed 's/Tumor_Seq_Allele2/ALT/' |
        sed 's/Tumor_Sample_Barcode/SAMPLE/'
        {params.additional_commands}
        > {output.maf}
        """)

ONCODRIVE_BUILD_DICT = {
    'hg38': 'hg38',
    'grch37': 'hg19'
}

def _get_region(wildcards):
    CFG = config["lcr-modules"]["oncodriveclustl"]
    build_regions = CFG["regions_files"][wildcards.genome_build]
    if wildcards.genome_build == "grch37":
        regions_file = reference_files(build_regions[wildcards.region])
    else:
        # hg38 regions file for Oncodrive are not available in their bitbucket and haven't been downloaded through reference files workflow
        regions_file = build_regions[wildcards.region]

    return regions_file

rule _oncodriveclustl_run:
    input:
        maf = str(rules._oncodriveclustl_format_input.output.maf),
        reference = lambda w: reference_files("downloads/oncodrive/{genome_build}/datasets/genomereference/" + ONCODRIVE_BUILD_DICT[w.genome_build] + ".master"),
        region = _get_region
    output:
        txt = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/elements_results.txt",
        tsv = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/clusters_results.tsv",
        png = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/quantile_quantile_plot.png"
    log:
        stdout = CFG["logs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodriveclustl.stdout.log",
        stderr = CFG["logs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodriveclustl.stderr.log"
    params:
        local_path = CFG["reference_files_directory"] + "{genome_build}/",
        build = lambda w: (w.genome_build).replace("grch37","hg19").replace("grch38","hg38"),
        command_line_options = CFG["options"]["clustl_options"] if CFG["options"]["clustl_options"] is not None else ""
    threads:
        CFG["threads"]["clustl"]
    resources:
        **CFG["resources"]["clustl"]
    conda: CFG["conda_envs"]["clustl"]
    shell:
        op.as_one_line("""
        export BGDATA_LOCAL={params.local_path} &&
        oncodriveclustl
        -i {input.maf}
        -o "$(dirname $(realpath {output.txt}))"
        -r {input.region}
        -g {params.build}
        --cores {threads}
        --qqplot
        {params.command_line_options}
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodriveclustl_out:
    input:
        txt = str(rules._oncodriveclustl_run.output.txt),
        tsv = str(rules._oncodriveclustl_run.output.tsv),
        png = str(rules._oncodriveclustl_run.output.png)
    output:
        txt = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/elements_results.txt",
        tsv = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/clusters_results.tsv",
        png = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/quantile_quantile_plot.png"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module=True)
        op.relative_symlink(input.txt, output.txt, in_module=True)
        op.relative_symlink(input.png, output.png, in_module=True)

rule _oncodriveclustl_get_cluster_coordinates:
    input:
        elements = str(rules._oncodriveclustl_run.output.txt),
        clusters = str(rules._oncodriveclustl_run.output.tsv)
    output:
        tsv = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/genomic_coordinates_clusters_results_{q_value}.tsv"
    params:
        q_value = lambda w: w.q_value,
        samples = CFG["detailed_clusters_options"]["minimum_samples"],
        p_value = CFG["detailed_clusters_options"]["p_value"],
        score = "--score " + CFG["detailed_clusters_options"]["minimum_score"] if CFG["detailed_clusters_options"]["minimum_score"] is not None else "",
        script = CFG["scripts"]["detailed_clusters_script"]
    log:
        stdout = CFG["logs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/detailed_clusters_{q_value}.stdout.log",
        stderr = CFG["logs"]["oncodriveclustl"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/detailed_clusters_{q_value}.stderr.log"
    shell:
        op.as_one_line("""
        {params.script}
        -e {input.elements}
        -c {input.clusters}
        -q {params.q_value}
        -n {params.samples}
        -p {params.p_value}
        -o {output.tsv}
        {params.score}
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodriveclustl_genomic_coordinates_out:
    input:
        genomic_coordinates = str(rules._oncodriveclustl_get_cluster_coordinates.output.tsv)
    output:
        genomic_coordinates = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/genomic_coordinates_clusters_results_{q_value}.tsv"
    run:
        op.relative_symlink(input.genomic_coordinates, output.genomic_coordinates)

def _get_oncodriveclustl_outputs(wildcards):
    CFG = config["lcr-modules"]["oncodriveclustl"]
    checkpoint_output = os.path.dirname(str(checkpoints._oncodriveclustl_prep_input.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output+"/{md5sum}.maf.content")

    return expand(
        [
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/elements_results.txt",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/clusters_results.tsv",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/quantile_quantile_plot.png",
        ],
        md5sum = SUMS
    )

def _get_oncodriveclustl_coordinates(wildcards):
    CFG = config["lcr-modules"]["oncodriveclustl"]
    checkpoint_output = os.path.dirname(str(checkpoints._oncodriveclustl_prep_input.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output+"/{md5sum}.maf.content")

    return expand(
        CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/genomic_coordinates_clusters_results_{q_value}.tsv",
        md5sum = SUMS,
        q_value = CFG["q_values"]
    )

rule _oncodriveclustl_aggregate:
    input:
        oncodrive_outputs = _get_oncodriveclustl_outputs,
        coordinate_output = _get_oncodriveclustl_coordinates
    output:
        aggregate = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/aggregate/{sample_set}--{launch_date}--{region}.done"
    shell:
        "touch {output.aggregate}"

# Generates the target sentinels for each run, which generate the symlinks
rule _oncodriveclustl_all:
    input:
        expand(
            [
                str(rules._oncodriveclustl_aggregate.output.aggregate)
            ],
            genome_build=CFG["genome_builds"],
            sample_set=CFG["maf_processing"]["sample_sets"],
            region = CFG["regions"],
            launch_date = launch_date
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
