#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "kallisto", 
    version = "1.0",
    subdirs = ["inputs", "kallisto_quant", "outputs"],
    req_references = ["kal_index", "kal_gtf", "kal_chrom"] 
)

localrules: 
    _kallisto_input,
    _kallisto_output,
    _kallisto_all


##### RULES #####
rule _kallisto_input:
    input:
        fq = CFG["inputs"]["fastq"]
    output:
        bfq = expand("{fqDIR}{{seq_type}}--{{genome_build}}/{{sample_id}}.{num}.fastq.gz", fqDIR = CFG["dirs"]["inputs"], num = ["1", "2"]) 
    run:
        md.symlink(input.fq[0], output.fq[0])
        md.symlink(input.fq[1], output.fq[1])


rule _kallisto_quant:
    input:
        fq = rules._kallisto_input.output.fq
    output:
        tsv = CFG["dirs"]["kallisto_quant"] + "{seq_type}--{genome_build}/{sample_id}/abundance.tsv"
    log:
        stdout = CFG["logs"]["kallisto_quant"] + "{seq_type}--{genome_build}/{sample_id}/quant.stdout.log",
        stderr = CFG["logs"]["kallisto_quant"] + "{seq_type}--{genome_build}/{sample_id}/quant.stderr.log"
    params:
        opts = CFG["options"]["quant"],
        strand = md.switch_on_column("strand", CFG["samples"], CFG["options"]["strand"]),
        idx = md.get_references(CFG, "kal_index"),
        outDIR = CFG["dirs"]["kallisto_quant"] + "{seq_type}--{genome_build}/{sample_id}/"
    conda:
        CFG["conda_envs"].get("kallisto") or "envs/kallisto-0.46.yaml"
    threads:
        CFG["threads"]["quant"]
    resources: 
        mem_mb = CFG["mem_mb"]["quant"]
    shell:
        md.as_one_line("""
        kallisto quant --threads {threads} 
        -i {params.idx} -o {params.outDIR}
        {params.strand}
        {params.opts}
        {input.fq}
        > {log.stdout} 2> {log.stderr}
        """)

rule _kallisto_output:
    input:
        tsv = rules._kallisto_quant.output.tsv
    output:
        tsv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.abundance.tsv"
    run:
        md.symlink(input.tsv, output.tsv)


rule _kallisto_all:
    input:
        expand(rules._kallisto_output.output.tsv, 
               zip,
               seq_type = CFG["samples"]["seq_type"],
               genome_build = CFG["samples"]["genome_build"],
               sample_id = CFG["samples"]["sample_id"])


##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
