#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["outputs"]`
CFG = op.setup_module(
    name = "bam2fastq",
    version = "1.0",
    subdirectories = ["inputs", "fastq", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _bam2fastq_input_bam,
    _bam2fastq_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _bam2fastq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)


rule _bam2fastq_run:
    input:
        bam = rules._bam2fastq_input_bam.output.bam
    output:
        fastq = expand("{fq_dir}{{seq_type}}--{{genome_build}}/{{sample_id}}.{read_num}.fastq", fq_dir = CFG["dirs"]["fastq"], read_num = ["read1", "read2"])
    log:
        stdout = CFG["logs"]["fastq"] + "{seq_type}--{genome_build}/{sample_id}/bam2fastq.stdout.log",
        stderr = CFG["logs"]["fastq"] + "{seq_type}--{genome_build}/{sample_id}/bam2fastq.stderr.log"
    params:
        opts = CFG["options"]["bam2fastq"]
    conda:
        CFG["conda_envs"]["picard"]
    threads:
        CFG["threads"]["bam2fastq"]
    resources:
        mem_mb = CFG["mem_mb"]["bam2fastq"]
    shell:
        op.as_one_line("""
        picard -Xmx{resources.mem_mb}m SamToFastq {params.opts}
        I={input.bam} FASTQ=>(gzip > {output.fq[0]}) SECOND_END_FASTQ=>(gzip > {output.fq[1]}) 
        > {log.stdout} &> {log.stderr}
        """)


rule _bam2fastq_output:
    input:
        fastq = rules._bam2fastq_run.output.fastq
    output:
        fastq_1 = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.read1.fastq",
        fastq_2 = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.read2.fastq"
       # expand("{fq_dir}{{seq_type}}--{{genome_build}}/{{sample_id}}.{read_num}.fastq", fq_dir = CFG["dirs"]["outputs"], read_num = ["read1", "read2"])
    run:
        op.relative_symlink(input.fastq[0], output.fastq_1)
        op.relative_symlink(input.fastq[1], output.fastq_2)

print(CFG["samples"]["sample_id"], type(CFG["samples"]["sample_id"]))

rule _bam2fastq_all:
    input:
        expand(
            [
                rules._bam2fastq_run.output.fastq,
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=["TCRBOA7-N-WEX", "TCRBOA7-T-WEX", "TCRBOA7-ALT-T-WEX", "TCRBOA7-T-RNA"])
            #CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
