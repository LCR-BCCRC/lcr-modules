#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Helena Winata
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["kallisto"]`
CFG = op.setup_module(
    name = "kallisto",
    version = "1.0",
    subdirectories = ["inputs", "kallisto", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _kallisto_input_fastq,
    _kallisto_output_tsv,
    _kallisto_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _kallisto_input_fastq:
    input:
        fastq_1 = CFG["inputs"]["sample_fastq_1"],
        fastq_2 = CFG["inputs"]["sample_fastq_2"]
    output:
        fastq = expand("{fqDIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{read_num}.fastq.gz", fqDIR = CFG["dirs"]["inputs"], read_num = ["read1", "read2"]) 
    run:
        op.relative_symlink(input.fastq_1, output.fastq[0])
        op.relative_symlink(input.fastq_2, output.fastq[1])


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _kallisto_quant:
    input:
        fastq = expand("{fqDIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{read_num}.fastq.gz", fqDIR = CFG["dirs"]["inputs"], read_num = ["read1", "read2"]) 
    output:
        tsv = CFG["dirs"]["kallisto"] + "{seq_type}--{genome_build}/{sample_id}/abundance.tsv"
    log:
        stdout = CFG["logs"]["kallisto"] + "{seq_type}--{genome_build}/{sample_id}/quant.stdout.log",
        stderr = CFG["logs"]["kallisto"] + "{seq_type}--{genome_build}/{sample_id}/quant.stderr.log"
    params:
        opts = CFG["options"]["quant"],
        strand = op.switch_on_column("strand", CFG["samples"], CFG["options"]["strand"], match_on="sample"),
        idx = reference_files("genomes/{genome_build}/kallisto_index/0.46.2/transcriptome.idx"),
        gtf = reference_files("genomes/{genome_build}/annotations/gencode_annotation-33.gtf"),
        chrom_sizes = reference_files("genomes/{genome_build}/genome_fasta/genome_chrom_sizes.txt")
    conda:
        CFG["conda_envs"]["kallisto"]
    threads:
        CFG["threads"]["quant"]
    resources:
        mem_mb = CFG["mem_mb"]["quant"]
    shell:
        op.as_one_line("""
        kallisto quant --threads {threads} 
        -i {params.idx} 
        -o $(dirname {output.tsv})
        {params.strand}
        {params.opts}
        {input.fastq}
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _kallisto_output_tsv:
    input:
        tsv = rules._kallisto_quant.output.tsv
    output:
        tsv = CFG["dirs"]["outputs"] + "counts/{seq_type}--{genome_build}/{sample_id}.abundance.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv)


# Generates the target sentinels for each run, which generate the symlinks
rule _kallisto_all:
    input:
        expand(
            [
                rules._kallisto_output_tsv.output.tsv,
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
