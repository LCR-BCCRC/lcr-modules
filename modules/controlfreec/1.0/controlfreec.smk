#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper
# Module Author:    Jasper
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`

# `CFG` is a shortcut to `config["lcr-modules"]["controlfreec"]`
CFG = op.setup_module(
    name = "controlfreec",
    version = "1.0",
    subdirectories = ["inputs", "run", "calc_sig", "plot", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _controlfreec_input_bam,
    _controlfreec_config,
    _controlfreec_run,
    _controlfreec_calc_sig,
    _controlfreec_plot,
    _controlfreec_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _controlfreec_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)

rule _controlfreec_config:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
    output:
        CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/config_WGS.txt"
    params:
        config = "config/freec/config_WGS.txt",
        outdir = CFG["dirs"]["run"] + "{sample_id}"
    shell:
        "sed \"s|BAMFILE|{input}|g\" {params.config} | "
        "sed \"s|OUTDIR|{params.outdir}|g\" > {output}"

rule _controlfreec_run:
    input:
        CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/config_WGS.txt"
    output:
        CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_info.txt",
        CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_ratio.txt",
        CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_CNVs"
    # params:
    # 	FREEC = config["software"]["FREEC"]
    conda: CFG["conda_envs"]["controlfreec"]
    threads: CFG["threads"]["controlfreec_run"]
    resources: mem_mb = CFG["mem_mb"]["controlfreec_run"]
    log:
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/run.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/run.stderr.log"
    shell:
        "freec -conf {input} > {log.stdout} 2> {log.stderr} "

rule _controlfreec_calc_sig:
    input:
        CNVs = CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_CNVs",
        ratios = CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_ratio.txt",
    output:
        CFG["dirs"]["calc_sig"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"],
        R = "/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/ichorcna/bin/R"
    threads: CFG["threads"]["calc_sig"]
    resources: mem_mb = CFG["mem_mb"]["calc_sig"]
    log:         
        stdout = CFG["logs"]["calc_sig"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stdout.log",
        stderr = CFG["logs"]["calc_sig"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stderr.log"
    shell:
        "cat {params.calc_sig} | {params.R} --slave --args {input.CNVs} {input.ratios} > {log.stdout} 2> {log.stderr}"

rule _controlfreec_plot:
    input:
        ratios = CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_ratio.txt",
        # BAF = CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_BAF.txt",
        info = CFG["dirs"]["run"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_info.txt"
    output:
        plot = CFG["dirs"]["plot"] + "{sample_id}/{seq_type}--{genome_build}/{sample_id}.bam_ratio.txt.png"
    params:
        plot = CFG["software"]["FREEC_graph"],
        R = "/projects/dscott_prj/CCSRI_1500/lowpassWGS/envs/ichorcna/bin/R"
    threads: CFG["threads"]["plot"]
    resources: mem_mb = CFG["mem_mb"]["plot"]
    log: 
        stdout = CFG["logs"]["plot"] + "{seq_type}--{genome_build}/{sample_id}/plot.stdout.log",
        stderr = CFG["logs"]["plot"] + "{seq_type}--{genome_build}/{sample_id}/plot.stderr.log"
    shell:
        "cat {params.plot} | {params.R} --slave --args `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} > {log.stdout} 2> {log.stderr} "


# Generates the target sentinels for each run, which generate the symlinks
rule _controlfreec_all:
    input:
        expand(
            [
                # str(rules._controlfreec_output_txt.output.txt),
                str(rules._controlfreec_plot.output.plot)
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
