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
    subdirectories = ["inputs", "run", "outputs"]
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
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/config_WGS.txt"
    conda:
        CFG["conda_envs"]["sambamba"]
    params:
        config = CFG["options"]["configFile"],
        outdir = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/",
        window = CFG["options"]["window"],
        ploidy = CFG["options"]["ploidy"],
        chrLen = CFG["options"]["chrLenFile"],
        chrFiles = CFG["options"]["chrFiles"],
        reference = CFG["options"]["gemMappabilityFile"]
    shell:
        "samtoolspath=$(which samtools ) ; "
        "samtoolsPathName=$(echo $samtoolspath) ; "
        "sambambapath=$(which sambamba ) ; "
        "sambambaPathName=$(echo $sambambapath) ; "
        "bedtoolspath=$(which bedtools ) ; "
        "bedtoolsPathName=$(echo $bedtoolspath) ; "
        "sed \"s|BAMFILE|{input.bam}|g\" {params.config} | "
        "sed \"s|OUTDIR|{params.outdir}|g\" | "
        "sed \"s|windowSize|{params.window}|g\" | "
        "sed \"s|ploidyInput|{params.ploidy}|g\" | "
        "sed \"s|chrLenFileInput|{params.chrLen}|g\" | "
        "sed \"s|chrFilesPath|{params.chrFiles}|g\" | "
        "sed \"s|sambambaPath|$sambambaPathName|g\" | "
        "sed \"s|bedtoolsPath|$bedtoolsPathName|g\" | "
        "sed \"s|samtoolsPath|$samtoolsPathName|g\" | "
        "sed \"s|referenceFile|{params.reference}|g\" > {output}"

rule _controlfreec_run:
    input:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/config_WGS.txt"
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_info.txt",
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs"
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
        CNVs = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs",
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
    output:
        CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_CNVs.p.value.txt"
    params:
        calc_sig = CFG["software"]["FREEC_sig"]
    threads: CFG["threads"]["calc_sig"]
    resources: mem_mb = CFG["mem_mb"]["calc_sig"]
    conda: CFG["conda_envs"]["controlfreec"]
    log:         
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/calc_sig.stderr.log"
    shell:
        "cat {params.calc_sig} | R --slave --args {input.CNVs} {input.ratios} > {log.stdout} 2> {log.stderr}"

rule _controlfreec_plot:
    input:
        ratios = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt",
        # BAF = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_BAF.txt",
        info = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_info.txt"
    output:
        plot = CFG["dirs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/{sample_id}.bam_ratio.txt.png"
    params:
        plot = CFG["software"]["FREEC_graph"]
    threads: CFG["threads"]["plot"]
    resources: mem_mb = CFG["mem_mb"]["plot"]
    conda: CFG["conda_envs"]["controlfreec"]
    log: 
        stdout = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/plot.stdout.log",
        stderr = CFG["logs"]["run"] + "{seq_type}--{genome_build}/{sample_id}/plot.stderr.log"
    shell:
        "cat {params.plot} | R --slave --args `grep \"Output_Ploidy\" {input.info} | cut -f 2` {input.ratios} > {log.stdout} 2> {log.stderr} "


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
