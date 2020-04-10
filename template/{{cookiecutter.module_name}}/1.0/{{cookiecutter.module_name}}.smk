#!/usr/bin/env snakemake


##### SETUP #####


# Import package with useful functions for developing analysis modules.
import modutils as md

# Make sure the `CFG` variable doesn't exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
# `CFG` corresponds to `config["lcr-modules"]["{{cookiecutter.module_name}}"]`
CFG = md.setup_module(
    config = config, 
    name = "{{cookiecutter.module_name}}", 
    version = "1.0",
    # TODO: If applicable, add other output subdirectories
    subdirs = ["inputs", "{{cookiecutter.module_name}}", "outputs"],
    # TODO: Replace "genome_fasta" with actual reference requirements
    req_references = ["genome_fasta"] 
)

# Define rules to be run locally when using a compute cluster.
# TODO: Replace with actual rule names
localrules: 
    _{{cookiecutter.module_name}}_input_bam,
    _{{cookiecutter.module_name}}_step_2,
    _{{cookiecutter.module_name}}_output_vcf,
    _{{cookiecutter.module_name}}_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/').
# TODO: Update input file (and if applicable, add one rule for each input file)
rule _{{cookiecutter.module_name}}_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)

{% if cookiecutter.module_run_per == "tumour" %}
# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_1:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"
    output:
        vcf = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
        fasta = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("{{cookiecutter.module_name}}") or "envs/{{cookiecutter.module_name}}.yaml"
    threads:
        CFG["threads"]["step_1"]
    resources: 
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        md.as_one_line("""
        <TODO> {params.opts} --tumour {input.tumour_bam} --normal {input.normal_bam} 
        --ref-fasta {params.fasta} --output {output.vcf} --threads {threads} 
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_2:
    input:
        vcf = rules._{{cookiecutter.module_name}}_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/').
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _{{cookiecutter.module_name}}_output_vcf:
    input:
        vcf = rules._{{cookiecutter.module_name}}_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.variants.filt.vcf"
    run:
        md.symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_{{cookiecutter.module_name}}_output_*` rule
rule _{{cookiecutter.module_name}}_all:
    input:
        expand(rules._{{cookiecutter.module_name}}_output_vcf.output.vcf, 
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
        bam = rules._{{cookiecutter.module_name}}_input_bam.output.bam
    output:
        vcf = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stdout.log",
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
        fasta = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("{{cookiecutter.module_name}}") or "envs/{{cookiecutter.module_name}}.yaml"
    threads:
        CFG["threads"]["step_1"]
    resources: 
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        md.as_one_line("""
        <TODO> {params.opts} --bam {input.bam} --ref-fasta {params.fasta} 
        --output {output.vcf} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _{{cookiecutter.module_name}}_step_2:
    input:
        vcf = rules._{{cookiecutter.module_name}}_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["{{cookiecutter.module_name}}"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/').
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _{{cookiecutter.module_name}}_output_vcf:
    input:
        vcf = rules._{{cookiecutter.module_name}}_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.variants.filt.vcf"
    run:
        md.symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_{{cookiecutter.module_name}}_output_*` rule
rule _{{cookiecutter.module_name}}_all:
    input:
        expand(rules._{{cookiecutter.module_name}}_output.output.vcf, 
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["samples"]["seq_type"],
               genome_build=CFG["samples"]["genome_build"],
               sample_id=CFG["samples"]["sample_id"])

{% endif %}
##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
