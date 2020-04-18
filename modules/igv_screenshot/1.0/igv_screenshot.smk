#!/usr/bin/env snakemake


##### SETUP #####


# Import package with useful functions for developing analysis modules.
import modutils as md

# Make sure the `CFG` variable doesn't exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
# `CFG` corresponds to `config["lcr-modules"]["igv_screenshot"]`
CFG = md.setup_module(
    config = config, 
    name = "igv_screenshot", 
    version = "1.0",
    # TODO: If applicable, add other output subdirectories
    subdirs = ["inputs", "igv_screenshot", "outputs"],
    # TODO: Replace "genome_fasta" with actual reference requirements
    req_references = ["genome_fasta"] 
)

# Define rules to be run locally when using a compute cluster.
# TODO: Replace with actual rule names
localrules: 
    _igv_screenshot_input_bam,
    _igv_screenshot_step_2,
    _igv_screenshot_output_vcf,
    _igv_screenshot_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/').
# TODO: Update input file (and if applicable, add one rule for each input file)
rule _igv_screenshot_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)

#
rule _igv_screenshot_input_maf:
    input:
        bam = CFG["inputs"]["sample_maf"]
    output:
        bam = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)

# generates batch files for IGV
rule generate_batch:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"
	    maf_file = MAF_DIR + "{id}/{id}.filt_common.maf",
        geneList = BED_DIR + "geneList.txt"
    output: 
        batchfile = BATCH_DIR + "{id}.igv.batch"
    params: 
        genome=GENOME_BUILD,
        pad=PADDING, 
        dirPic= SNAPSHOT_DIR + "{id}/", 
        mode="single"
    conda: "igvTest.yml"
    threads: 1
    resources:
        mem_mb = 8000
    message:
        "Generating batch files....."
    shell:
        #'ls {input.maf_file} {input.nbam} {input.nbai} {wildcards.id} > {output.batchfile} && touch {output.batchfile}'
        'python3 generate_igv_batch4.py --genome_build {params.genome} --input {input.maf_file} --mode {params.mode} --bam {input.nbam} --basename {wildcards.id} --snapshot_dir {params.dirPic} --pad {params.pad} --batchfile {output.batchfile} --geneList {input.geneList}'


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _igv_screenshot_step_1:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam"
    output:
        vcf = CFG["dirs"]["igv_screenshot"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.vcf.gz"
    log:
        stdout = CFG["logs"]["igv_screenshot"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
        stderr = CFG["logs"]["igv_screenshot"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
    params:
        opts = CFG["options"]["step_1"]
        fasta = lambda wildcards: config["reference"][wildcards.genome_build]["genome_fasta"]
    conda:
        CFG["conda_envs"].get("igv_screenshot") or "envs/igv_screenshot.yaml"
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
rule _igv_screenshot_step_2:
    input:
        vcf = rules._igv_screenshot_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["igv_screenshot"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["igv_screenshot"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/').
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _igv_screenshot_output_vcf:
    input:
        vcf = rules._igv_screenshot_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.variants.filt.vcf"
    run:
        md.symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_igv_screenshot_output_*` rule
rule _igv_screenshot_all:
    input:
        expand(rules._igv_screenshot_output_vcf.output.vcf, 
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
