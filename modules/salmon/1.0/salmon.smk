#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Helena Winata
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["salmon"]`
CFG = op.setup_module(
    name = "salmon",
    version = "1.0",
    subdirectories = ["inputs", "salmon", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _salmon_input_fastq,
    _salmon_output,
    _salmon_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _salmon_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"]
    output:
        fastq_1 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R1.fastq",
        fastq_2 = CFG["dirs"]["inputs"] + "fastq/{seq_type}--{genome_build}/{sample_id}.R2.fastq"
    run:
        op.relative_symlink(input.fastq_1, output.fastq_1)
        op.relative_symlink(input.fastq_2, output.fastq_2)


rule _salmon_quant:
    input:
        fastq_1 = rules._salmon_input_fastq.output.fastq_1,
        fastq_2 = rules._salmon_input_fastq.output.fastq_2,
        index = reference_files("genomes/{genome_build}/salmon_index/salmon-1.3.0/index")
    output:
        quant = CFG["dirs"]["salmon"] + "{seq_type}--{genome_build}/{sample_id}/quant.sf"
    log:
        stdout = CFG["logs"]["salmon"] + "{seq_type}--{genome_build}/{sample_id}/quant.stdout.log",
        stderr = CFG["logs"]["salmon"] + "{seq_type}--{genome_build}/{sample_id}/quant.stderr.log"
    params:
        opts = CFG["options"]["quant"]
    conda:
        CFG["conda_envs"]["salmon"]
    threads:
        CFG["threads"]["quant"]
    resources:
        mem_mb = CFG["mem_mb"]["quant"]
    shell:
        op.as_one_line("""
        salmon quant -p {threads} 
        {params.opts} 
        -i {input.index} 
        -o $(dirname {output.quant})
        -1 {input.fastq_1} -2 {input.fastq_2} 
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _salmon_output:
    input:
        quant = rules._salmon_quant.output.quant
    output:
        quant = CFG["dirs"]["outputs"] + "quant/{seq_type}--{genome_build}/{sample_id}.quant.sf"
    run:
        op.relative_symlink(input.quant, output.quant)


# Generates the target sentinels for each run, which generate the symlinks
rule _salmon_all:
    input:
        expand(rules._salmon_output.output.quant,
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
