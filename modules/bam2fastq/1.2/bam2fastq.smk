#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["bam2fastq"]`
CFG = op.setup_module(
    name = "bam2fastq",
    version = "1.2",
    subdirectories = ["inputs", "fastq", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _bam2fastq_input_bam,
    _bam2fastq_output,
    _bam2fastq_all,

assert type(CFG["temp_outputs"])==bool, (
    "The config value for temp_outputs must be set to a boolean value. "
    "Note that TRUE, FALSE, \"True\", and \"False\" are not Python booleans. "
)

##### RULES #####

# Functions to obtain the genome_build per sample from the sample table
def get_genome_build(wildcards):
    tbl = config["lcr-modules"]["bam2fastq"]["samples"]
    genome_build = tbl[(tbl.sample_id == wildcards.sample_id) & (tbl.seq_type == wildcards.seq_type)]["genome_build"]
    return reference_files(expand("genomes/{genome_build}/genome_fasta/genome.fa", genome_build = genome_build))

def get_bams(wildcards,build = False):
    CFG = config["lcr-modules"]["bam2fastq"]
    tbl = CFG["samples"]
    return(expand(
        CFG["dirs"]["inputs"] + "{{seq_type}}--{genome_build}/{{sample_id}}.bam", 
        genome_build = tbl[(tbl.sample_id == wildcards.sample_id) & (tbl.seq_type == wildcards.seq_type)]["genome_build"]
    ))


rule _bam2fastq_input_bam:
    input:
        bam_path = CFG['inputs']['sample_bam']
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input, output.bam)


# Conditional rules depending on whether or not fastq outputs will be temporary
if CFG["temp_outputs"]:
    rule _bam2fastq_run:
        input:
            bam = get_bams,
            genome = get_genome_build
        output:
            fastq_1 = temp(CFG["dirs"]["fastq"] + "{seq_type}/{sample_id}.read1.fastq.gz"),
            fastq_2 = temp(CFG["dirs"]["fastq"] + "{seq_type}/{sample_id}.read2.fastq.gz")
        log:
            stdout = CFG["logs"]["fastq"] + "{seq_type}/{sample_id}/bam2fastq.stdout.log",
            stderr = CFG["logs"]["fastq"] + "{seq_type}/{sample_id}/bam2fastq.stderr.log"
        params:
            opts = CFG["options"]["bam2fastq"]
        conda:
            CFG["conda_envs"]["picard"]
        threads:
            CFG["threads"]["bam2fastq"]
        resources:
            **CFG["resources"]["bam2fastq"]
        shell:
            op.as_one_line("""
            picard -Xmx{resources.mem_mb}m SamToFastq {params.opts}
            I={input.bam} FASTQ=>(gzip > {output.fastq_1}) SECOND_END_FASTQ=>(gzip > {output.fastq_2}) 
            REFERENCE_SEQUENCE={input.genome}
            > {log.stdout} &> {log.stderr}
            """)

elif not CFG["temp_outputs"]:
    rule _bam2fastq_run:
        input:
            bam = get_bams,
            genome = get_genome_build
        output:
            fastq_1 = CFG["dirs"]["fastq"] + "{seq_type}/{sample_id}.read1.fastq.gz",
            fastq_2 = CFG["dirs"]["fastq"] + "{seq_type}/{sample_id}.read2.fastq.gz"
        log:
            stdout = CFG["logs"]["fastq"] + "{seq_type}/{sample_id}/bam2fastq.stdout.log",
            stderr = CFG["logs"]["fastq"] + "{seq_type}/{sample_id}/bam2fastq.stderr.log"
        params:
            opts = CFG["options"]["bam2fastq"]
        conda:
            CFG["conda_envs"]["picard"]
        threads:
            CFG["threads"]["bam2fastq"]
        resources:
            **CFG["resources"]["bam2fastq"]
        shell:
            op.as_one_line("""
            picard -Xmx{resources.mem_mb}m SamToFastq {params.opts}
            I={input.bam} FASTQ=>(gzip > {output.fastq_1}) SECOND_END_FASTQ=>(gzip > {output.fastq_2}) 
            REFERENCE_SEQUENCE={input.genome}
            > {log.stdout} &> {log.stderr}
            """)
else:
    raise ValueError("CFG['temp_outputs'] must be set to a boolean value (True or False)")


rule _bam2fastq_output:
    input:
        fastq_1 = str(rules._bam2fastq_run.output.fastq_1),
        fastq_2 = str(rules._bam2fastq_run.output.fastq_2)
    output:
        fastq_1 = CFG["dirs"]["outputs"] + "{seq_type}/{sample_id}.read1.fastq.gz",
        fastq_2 = CFG["dirs"]["outputs"] + "{seq_type}/{sample_id}.read2.fastq.gz"
    run:
        op.relative_symlink(input.fastq_1, output.fastq_1)
        op.relative_symlink(input.fastq_2, output.fastq_2)


rule _bam2fastq_all:
    input:
        expand(
            [
                rules._bam2fastq_output.output.fastq_1,
                rules._bam2fastq_output.output.fastq_2,
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            sample_id=CFG["samples"]["sample_id"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
