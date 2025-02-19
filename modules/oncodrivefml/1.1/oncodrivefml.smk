#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Manuela Cruz
# Module Author:    Manuela Cruz
# Contributors:     N/A

##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
from datetime import datetime
import numpy as np

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["oncodrivefml"]`
CFG = op.setup_module(
    name = "oncodrivefml",
    version = "1.1",
    subdirectories = ["inputs", "prepare_mafs", "cadd", "oncodrivefml", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _oncodrivefml_input_maf,
    _oncodrivefml_sample_set,
    _oncodrivefml_prep_input,
    _oncodrivefml_blacklist,
    _oncodrivefml_format_input,
    _oncodrivefml_get_hg19_scores,
    _oncodrivefml_unzip_output,
    _oncodrivefml_out,
    _oncodrivefml_symlink_content,
    _oncodrivefml_aggregate,
    _oncodrivefml_all,


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG["launch_date"]
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later
PREPARE_MAFS = os.path.abspath(config["lcr-modules"]["oncodrivefml"]["maf_processing"]["prepare_mafs"])

# Symlink master MAF to input directory, easier to extract files this way
rule _oncodrivefml_input_maf:
    input:
        maf = CFG["inputs"]["input_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "master_maf/{seq_type}--{genome_build}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

rule _oncodrivefml_sample_set:
    input:
        subsetting_categories = ancient(CFG["inputs"]["subsetting_categories"])
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

checkpoint _oncodrivefml_prep_input:
    input:
        maf = expand(
            str(rules._oncodrivefml_input_maf.output.maf),
            seq_type = CFG["samples"]["seq_type"].unique(),
            allow_missing=True),
        subsetting_categories = ancient(str(rules._oncodrivefml_sample_set.output.subsetting_categories))
    output:
        CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/done"
    log:
        stdout = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/prep_input/prep_input_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["maf_processing"]["include_non_coding"]).upper(),
        mode = "OncodriveFML",
        metadata_cols = CFG["samples"],
        metadata_dim = CFG["samples"].shape,
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

rule _oncodrivefml_blacklist:
    input:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.maf",
        blacklists = CFG["maf_processing"]["blacklists"],
        deblacklist_script = CFG["scripts"]["deblacklist_script"]
    output:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.deblacklisted.maf"
    log:
        stdout = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stdout.log",
        stderr = CFG["logs"]["prepare_mafs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stderr.log"
    shell:
        op.as_one_line("""
        {input.deblacklist_script}
        --input {input.maf}
        --output {output.maf}
        --blacklists {input.blacklists}
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodrivefml_format_input:
    input:
        maf = str(rules._oncodrivefml_blacklist.output.maf)
    output:
        maf = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.fml_input.maf"
    params:
        columns = CFG["format_fml_input"]["input_columns"],
        additional_commands = CFG["format_fml_input"]["additional_commands"]
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
    'grch38': 'hg38',
    'grch37': 'hg19',
    'hg38': 'hg38',
    'hg19': 'hg19'
}

REFERENCE_FILES_VERSION_DICT = {
    'hg38': 'grch38',
    'grch38': 'grch38',
    'grch37': 'grch37',
    'hg19': 'grch37'
}

rule _oncodrivefml_get_hg19_scores:
    output:
        scores = CFG["dirs"]["cadd"] + "grch37/genomicscores/caddpack/1.0-20170217/whole_genome_SNVs.fml",
        index = CFG["dirs"]["cadd"] + "grch37/genomicscores/caddpack/1.0-20170217/whole_genome_SNVs.fml.idx"
    params:
        score_dir = CFG["dirs"]["cadd"] + "grch37/"
    conda: CFG["conda_envs"]["fml"]
    shell:
        op.as_one_line("""
        export BGDATA_LOCAL={params.score_dir} &&
        bgdata get genomicscores/caddpack/1.0
        """)

def _get_score_path(wildcards):
    # Use score path provided in config if available, otherwise use score provided by bbglab (grch37 only)
    CFG = config["lcr-modules"]["oncodrivefml"]
    if CFG["options"]["score_path"] is not None:
        score_path = CFG["options"]["score_path"]
    elif ONCODRIVE_BUILD_DICT[wildcards.genome_build] == "hg19":
        score_path = str(rules._oncodrivefml_get_hg19_scores.output.scores)
    return score_path

def _get_region(wildcards):
    CFG = config["lcr-modules"]["oncodrivefml"]
    build_regions = CFG["regions_files"][REFERENCE_FILES_VERSION_DICT[wildcards.genome_build]]
    if wildcards.genome_build in ["grch37","hg19"] and wildcards.region != "custom":
        regions_file = reference_files(build_regions[wildcards.region])
    else:
        # hg38 regions file / custom regions files for Oncodrive are not available in their bitbucket and haven't been downloaded through reference files workflow
        regions_file = build_regions[wildcards.region]
    return regions_file

rule _oncodrivefml_run:
    input:
        maf = str(rules._oncodrivefml_format_input.output.maf),
        config = CFG["options"]["config_path"],
        reference = lambda w: reference_files("downloads/oncodrive/" + REFERENCE_FILES_VERSION_DICT[w.genome_build] + "/datasets/genomereference/" + ONCODRIVE_BUILD_DICT[w.genome_build] + ".master"),
        scores = _get_score_path,
        region = _get_region
    output:
        tsv = CFG["dirs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/{md5sum}-oncodrivefml.tsv.gz",
        png = CFG["dirs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/{md5sum}-oncodrivefml.png",
        html = CFG["dirs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/{md5sum}-oncodrivefml.html"
    log:
        stdout = CFG["logs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.stdout.log",
        stderr = CFG["logs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.stderr.log"
    params:
        local_path = lambda w: config["lcr-modules"]["oncodrivefml"]["reference_files_directory"] + REFERENCE_FILES_VERSION_DICT[w.genome_build] + "/",
        command_line_options = CFG["options"]["fml_options"] if CFG["options"]["fml_options"] is not None else ""
    threads:
        CFG["threads"]["fml"]
    resources:
        **CFG["resources"]["fml"]
    conda: CFG["conda_envs"]["fml"]
    shell:
        op.as_one_line("""
        export BGDATA_LOCAL={params.local_path} &&
        export BGDATA_OFFLINE=TRUE &&
        oncodrivefml
        -i {input.maf}
        -o "$(dirname $(realpath {output.tsv}))"
        -e {input.region}
        -c {input.config}
        --cores {threads}
        {params.command_line_options}
        > {log.stdout} 2> {log.stderr}
        """)

rule _oncodrivefml_unzip_output:
    input:
        tsv = str(rules._oncodrivefml_run.output.tsv)
    output:
        tsv = CFG["dirs"]["oncodrivefml"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/{md5sum}--oncodrivefml.tsv"
    shell:
        "zcat {input.tsv} > {output.tsv}"

rule _oncodrivefml_out:
    input:
        gz = str(rules._oncodrivefml_run.output.tsv),
        png = str(rules._oncodrivefml_run.output.png),
        html = str(rules._oncodrivefml_run.output.html),
        tsv = str(rules._oncodrivefml_unzip_output.output.tsv)
    output:
        gz = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.tsv.gz",
        png = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.png",
        html = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.html",
        tsv = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/oncodrivefml.tsv"
    run:
        op.relative_symlink(input.gz, output.gz, in_module=True)
        op.relative_symlink(input.png, output.png, in_module=True)
        op.relative_symlink(input.html, output.html, in_module=True)
        op.relative_symlink(input.tsv, output.tsv, in_module=True)

rule _oncodrivefml_symlink_content:
    input:
        content = CFG["dirs"]["prepare_mafs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.maf.content",
        tsv = str(rules._oncodrivefml_run.output.tsv),
        png = str(rules._oncodrivefml_run.output.png),
        html = str(rules._oncodrivefml_run.output.html)
    output:
        content = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{region}/{md5sum}.maf.content"
    run:
        op.relative_symlink(input.content, output.content, in_module=True)

def _get_oncodrivefml_outputs(wildcards):
    CFG = config["lcr-modules"]["oncodrivefml"]
    checkpoint_output = os.path.dirname(str(checkpoints._oncodrivefml_prep_input.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output+"/{md5sum}.maf.content")

    return expand(
        [
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/oncodrivefml.tsv.gz",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/oncodrivefml.png",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/oncodrivefml.html",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/oncodrivefml.tsv",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/{{region}}/{md5sum}.maf.content"
        ],
        md5sum = SUMS
    )

rule _oncodrivefml_aggregate:
    input:
        _get_oncodrivefml_outputs
    output:
        aggregate = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/aggregate/{sample_set}--{launch_date}--{region}.done"
    shell:
        "touch {output.aggregate}"

rule _oncodrivefml_all:
    input:
        expand(
            str(rules._oncodrivefml_aggregate.output.aggregate),
            genome_build = CFG["genome_builds"],
            sample_set = CFG["maf_processing"]["sample_sets"],
            region = CFG["regions"],
            launch_date = launch_date
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
