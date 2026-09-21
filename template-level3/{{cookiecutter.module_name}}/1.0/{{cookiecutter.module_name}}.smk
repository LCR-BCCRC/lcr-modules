#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  {{cookiecutter.original_author}}
# Module Author:    {{cookiecutter.module_author}}
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import sys, os
from os.path import join
import glob
import hashlib
import oncopipe as op
from datetime import datetime
import numpy as np

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
# `CFG` is a shortcut to `config["lcr-modules"]["{{cookiecutter.module_name}}"]`
CFG = op.setup_module(
    name = "{{cookiecutter.module_name}}",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "prepare_{{cookiecutter.input_file_type}}", "{{cookiecutter.module_name}}", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}},
    _{{cookiecutter.module_name}}_input_subsetting_categories,
    _{{cookiecutter.module_name}}_prepare_{{cookiecutter.input_file_type}},
    _{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}},
    _{{cookiecutter.module_name}}_all,


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to the prepare script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_SCRIPT =  os.path.abspath(config["lcr-modules"]["{{cookiecutter.module_name}}"]["prepare_script"])

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}}:
    input:
        {{cookiecutter.input_file_type}} = CFG["inputs"]["master_{{cookiecutter.input_file_type}}"]
    output:
        {{cookiecutter.input_file_type}} = CFG["dirs"]["inputs"] + "{{cookiecutter.input_file_type}}/{seq_type}/{sample_set}--{projection}--{launch_date}/input.{{cookiecutter.input_file_type}}"
    run:
        op.absolute_symlink(input.{{cookiecutter.input_file_type}}, output.{{cookiecutter.input_file_type}})

# Symlinks the subsetting categories input file into the module results directory (under '00-inputs/')
rule _{{cookiecutter.module_name}}_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

# Prepare the {{cookiecutter.input_file_type}} file
checkpoint _{{cookiecutter.module_name}}_prepare_{{cookiecutter.input_file_type}}:
    input:
        {{cookiecutter.input_file_type}} = expand(
                    str(rules._{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}}.output.{{cookiecutter.input_file_type}}),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._{{cookiecutter.module_name}}_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["prepare_{{cookiecutter.input_file_type}}"] + "{sample_set}--{projection}--{launch_date}/done"
    log:
        CFG["logs"]["prepare_{{cookiecutter.input_file_type}}"] + "{sample_set}--{projection}--{launch_date}/prepare_{{cookiecutter.input_file_type}}.log"
    conda:
        CFG["conda_envs"]["prepare"]
    params:
        {% if cookiecutter.input_file_type == "maf" %}
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        {% endif %}
        # TODO: remove blank lines here (left here by the cookiecutter template)
        mode = "<TODO> must match what you used in the generate_smg_inputs lcr-script",
        metadata_cols = CFG["samples"],
        metadata_dim = CFG["samples"].shape,
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_SCRIPT


# Run {{cookiecutter.module_name}}
# TODO: Add all expected output files below
rule _{{cookiecutter.module_name}}_run:
    input:
        {{cookiecutter.input_file_type}} = CFG["dirs"]["prepare_{{cookiecutter.input_file_type}}"] + "{sample_set}--{projection}--{launch_date}/{md5sum}.{{cookiecutter.input_file_type}}"
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{sample_set}--{projection}/{launch_date}--{md5sum}/<TODO>.{{cookiecutter.output_file_type}}"
    log:
        stdout = CFG["logs"]["{{cookiecutter.module_name}}"] + "{sample_set}--{projection}/{launch_date}--{md5sum}/{{cookiecutter.module_name}}.stdout.log",
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{sample_set}--{projection}/{launch_date}--{md5sum}/{{cookiecutter.module_name}}.stderr.log"
    params:
        opts = CFG["options"]["{{cookiecutter.module_name}}_run"]
    conda:
        CFG["conda_envs"]["{{cookiecutter.module_name}}"]
    threads:
        CFG["threads"]["{{cookiecutter.module_name}}_run"]
    resources:
        **CFG["resources"]["{{cookiecutter.module_name}}_run"]
    shell:
        op.as_one_line("""
         <TODO> {params.opts} --input {input.{{cookiecutter.input_file_type}}} --output {output.{{cookiecutter.output_file_type}}} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add lines for each file meant to be exposed to the user
rule _{{cookiecutter.module_name}}_output:
    input:
        {{cookiecutter.output_file_type}} = str(rules._{{cookiecutter.module_name}}_run.output.{{cookiecutter.output_file_type}})
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["outputs"] + "{sample_set}--{projection}/{launch_date}--{md5sum}/<TODO>.{{cookiecutter.output_file_type}}"
    run:
        op.relative_symlink(input.{{cookiecutter.output_file_type}}, output.{{cookiecutter.output_file_type}}, in_module= True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["{{cookiecutter.module_name}}"]
    checkpoint_output = os.path.dirname(str(checkpoints._{{cookiecutter.module_name}}_prepare_{{cookiecutter.input_file_type}}.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum}.{{cookiecutter.input_file_type}}")
    return expand(
        [
            CFG["dirs"]["outputs"] + "{{ "{{" }}sample_set{{ "}}" }}--{{ "{{" }}projection{{ "}}" }}/{{ "{{" }}launch_date{{ "}}" }}--{md5sum}/<TODO>.{{cookiecutter.output_file_type}}",

        ],
        md5sum = SUMS
        )

# Aggregates outputs to remove md5sum from rule all
rule _{{cookiecutter.module_name}}_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{sample_set}--{projection}--{launch_date}.done"
    shell:
        op.as_one_line("""touch {output.aggregate}""")


rule _{{cookiecutter.module_name}}_all:
    input:
        expand(
            [
                CFG["dirs"]["prepare_{{cookiecutter.input_file_type}}"] + "{sample_set}--{projection}--{launch_date}/done",
                str(rules._{{cookiecutter.module_name}}_aggregate.output.aggregate)
            ],
            projection = CFG["projections"],
            sample_set = CFG["sample_set"],
            launch_date = launch_date
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)

