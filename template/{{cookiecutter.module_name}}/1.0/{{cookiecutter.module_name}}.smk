#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  {{cookiecutter.original_author}}
# Module Author:    {{cookiecutter.module_author}}
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["{{cookiecutter.module_name}}"]`
CFG = op.setup_module(
    name = "{{cookiecutter.module_name}}",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "{{cookiecutter.module_name}}", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}},
    _{{cookiecutter.module_name}}_step_2,
    _{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}},
    _{{cookiecutter.module_name}}_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
rule _{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}}:
    input:
        {{cookiecutter.input_file_type}} = CFG["inputs"]["sample_{{cookiecutter.input_file_type}}"]
    output:
        {{cookiecutter.input_file_type}} = CFG["dirs"]["inputs"] + "{{cookiecutter.input_file_type}}/{seq_type}--{genome_build}/{sample_id}.{{cookiecutter.input_file_type}}"
    run:
        op.relative_symlink(input.{{cookiecutter.input_file_type}}, output.{{cookiecutter.input_file_type}})

{% if cookiecutter.module_run_per == "tumour" %}
# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_1:
    input:
        tumour_{{cookiecutter.input_file_type}} = CFG["dirs"]["inputs"] + "{{cookiecutter.input_file_type}}/{seq_type}--{genome_build}/{tumour_id}.{{cookiecutter.input_file_type}}",
        normal_{{cookiecutter.input_file_type}} = CFG["dirs"]["inputs"] + "{{cookiecutter.input_file_type}}/{seq_type}--{genome_build}/{normal_id}.{{cookiecutter.input_file_type}}",
        fasta = reference_files(CFG["reference"]["genome_fasta"])
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.{{cookiecutter.output_file_type}}"
    log:
        stdout = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --tumour {input.tumour_{{cookiecutter.input_file_type}}} --normal {input.normal_{{cookiecutter.input_file_type}}}
        --ref-fasta {params.fasta} --output {output.{{cookiecutter.output_file_type}}} --threads {threads}
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_2:
    input:
        {{cookiecutter.output_file_type}} = rules._{{cookiecutter.module_name}}_step_1.output.{{cookiecutter.output_file_type}}
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.{{cookiecutter.output_file_type}}"
    log:
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.{{cookiecutter.output_file_type}}} > {output.{{cookiecutter.output_file_type}}} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}}:
    input:
        {{cookiecutter.output_file_type}} = rules._{{cookiecutter.module_name}}_step_2.output.{{cookiecutter.output_file_type}}
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["outputs"] + "{{cookiecutter.output_file_type}}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.{{cookiecutter.output_file_type}}"
    run:
        op.relative_symlink(input.{{cookiecutter.output_file_type}}, output.{{cookiecutter.output_file_type}})


# Generates the target sentinels for each run, which generate the symlinks
rule _{{cookiecutter.module_name}}_all:
    input:
        expand(
            [
                rules._{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}}.output.{{cookiecutter.output_file_type}},
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])

{% elif cookiecutter.module_run_per == "sample" %}
# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_1:
    input:
        {{cookiecutter.input_file_type}} = rules._{{cookiecutter.module_name}}_input_{{cookiecutter.input_file_type}}.output.{{cookiecutter.input_file_type}},
        fasta = reference_files(CFG["reference"]["genome_fasta"])
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/output.{{cookiecutter.output_file_type}}"
    log:
        stdout = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["step_1"]
    resources:
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        op.as_one_line("""
        <TODO> {params.opts} --input {input.{{cookiecutter.input_file_type}}} --ref-fasta {params.fasta}
        --output {output.{{cookiecutter.output_file_type}}} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_2:
    input:
        {{cookiecutter.output_file_type}} = rules._{{cookiecutter.module_name}}_step_1.output.{{cookiecutter.output_file_type}}
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/output.filt.{{cookiecutter.output_file_type}}"
    log:
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.{{cookiecutter.output_file_type}}} > {output.{{cookiecutter.output_file_type}}} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}}:
    input:
        {{cookiecutter.output_file_type}} = rules._{{cookiecutter.module_name}}_step_2.output.{{cookiecutter.output_file_type}}
    output:
        {{cookiecutter.output_file_type}} = CFG["dirs"]["outputs"] + "{{cookiecutter.output_file_type}}/{seq_type}--{genome_build}/{sample_id}.output.filt.{{cookiecutter.output_file_type}}"
    run:
        op.relative_symlink(input.{{cookiecutter.output_file_type}}, output.{{cookiecutter.output_file_type}})


# Generates the target sentinels for each run, which generate the symlinks
rule _{{cookiecutter.module_name}}_all:
    input:
        expand(
            [
                rules._{{cookiecutter.module_name}}_output_{{cookiecutter.output_file_type}}.output.{{cookiecutter.output_file_type}},
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])

{% endif %}
##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
