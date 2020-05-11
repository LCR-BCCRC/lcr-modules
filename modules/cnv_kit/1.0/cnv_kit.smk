#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

# Make sure the `CFG` variable doesn't exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "cnvkit", 
    version = "1.0",
    subdirs = ["inputs", "cnvkit", "outputs"],
    req_references = ["genome_fasta"] 
)

localrules: 
    _cnvkit_input,
    _cnvkit_output,
    _cnvkit_all


##### RULES #####

rule _cnvkit_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _cnvkit_batch:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam"
    output:
        cnr = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
        fasta = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("cnvkit") or "envs/cnvkit.yaml"
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
rule _cnvkit_step_2:
    input:
        vcf = rules._cnvkit_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["cnvkit"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/').
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _cnvkit_output_vcf:
    input:
        vcf = rules._cnvkit_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.variants.filt.vcf"
    run:
        md.symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_cnvkit_output_*` rule
rule _cnvkit_all:
    input:
        expand(rules._cnvkit_output_vcf.output.vcf, 
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["runs"]["tumour_seq_type"],
               genome_build=CFG["runs"]["tumour_genome_build"],
               tumour_id=CFG["runs"]["tumour_sample_id"],
               normal_id=CFG["runs"]["normal_sample_id"],
               pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
