#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "mutect", 
    version = "1.0",
    subdirs = ["inputs", "validate", "mutect", "filter", "outputs"],
    req_references = ["genome_fasta"] 
)

localrules: 
    _mutect_input,
    _mutect_output,
    _mutect_all


##### RULES #####

rule _mutect_input:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        md.symlink(input.bam, output.bam)
        md.symlink(input.bam + ".bai", output.bam + ".bai")


rule _mutect_validate:
    input:
        bam = rules._mutect_input.output.bam
    output:
        summary = CFG["dirs"]["validate"] + "{seq_type}--{genome_build}/{sample_id}.validate.txt"
    log:
        CFG["logs"]["validate"] + "{seq_type}--{genome_build}/{sample_id}/validate.stderr.log"
    params:
        opts = CFG["options"]["validate"],
        fasta = md.get_reference(CFG, "genome_fasta")
    conda:
        CFG["conda_envs"].get("gatk") or "envs/gatk-4.1.7.0.yaml"
    threads:
        CFG["threads"].get("validate") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("validate") or 4000
    shell:
        md.as_one_line("""
        gatk ValidateSamFile
        -I {input.bam}
        {params.opts}
        -O {output.summary}
        -R {params.fasta}
        2> {log}
        """)

'''
def _get_opt_output(wildcards, output):
    files = md.switch_on_wildcard("pair_status", {"_default": "--germline-resource {dir}/af-only-gnomad.vcf.gz --panel-of-normals {dir}/pon.vcf.gz", "no_normal": ""})
    print(os.path.dirname(output.vcf))
    paths = expand(files, dir = os.path.dirname(output.vcf))
    return paths
'''

rule _mutect_run:
    input:
        unpack(md.switch_on_wildcard("pair_status", {
            "_default" : {"normal_bam": CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam", "normal_sum" : CFG["dirs"]["validate"] + "{seq_type}--{genome_build}/{normal_id}.validate.txt"}, 
            "no_normal" : {}
            })),
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        tumour_summary = CFG["dirs"]["validate"] + "{seq_type}--{genome_build}/{tumour_id}.validate.txt"
    output:
        vcf = CFG["dirs"]["mutect"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.vcf.gz"
        #md.switch_on_wildcard("pair_status", {"_default" : CFG["dirs"]["mutect"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic.vcf.gz", "no_normal" : CFG["dirs"]["mutect"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/tumour.vcf.gz"})
    log:
        stdout = CFG["logs"]["mutect"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run.stdout.log",
        stderr = CFG["logs"]["mutect"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/run.stderr.log"
    params:
        opts = "", #md.switch_on_wildcard("pair_status", CFG["options"]["mutect"]),
        #_get_opt_output,
        #
        fasta = md.get_reference(CFG, "genome_fasta"),
        normal = md.switch_on_wildcard("pair_status", {
            "_default" : "-I {input.normal_bam} -normal {normal_id}",
            "no_normal" : ""}),
    conda:
        CFG["conda_envs"].get("gatk") or "envs/gatk-4.1.7.0.yaml"
    threads:
        CFG["threads"].get("mutect") or 4
    resources: 
        mem_mb = CFG["mem_mb"].get("mutect") or 6000
    shell:
        md.as_one_line("""
        gatk Mutect2 -R {params.fasta}
        -I {input.tumour_bam} 
        {params.normal}
        {params.opts}
        -O {output.vcf}
        > {log.stdout} 2> {log.stderr}
        """)


rule _mutect_filter:
    input:
        vcf = rules._mutect_run.output.vcf
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/variants.filtered.vcf.gz"
    log:
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter.stdout.log",
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter.stderr.log"
    params:
        opts = CFG["options"]["filter"],
        fasta = md.get_reference(CFG, "genome_fasta"),
    conda:
        CFG["conda_envs"].get("gatk") or "envs/gatk-4.1.7.0.yaml"
    threads:
        CFG["threads"].get("filter") or 2
    resources: 
        mem_mb = CFG["mem_mb"].get("filter") or 4000
    shell:
        md.as_one_line("""
        gatk FilterMutectCalls -R {params.fasta}
        -V {input.vcf} 
        {params.opts}
        -O {output.vcf}
        > {log.stdout} 2> {log.stderr}
        """)


rule _mutect_output:
    input:
        vcf = rules._mutect_filter.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.variants.filtered.vcf.gz"
    run:
        md.symlink(input.vcf, output.vcf)


rule _mutect_all:
    input:
        expand(rules._mutect_output.output.vcf, 
               zip,
               seq_type=CFG["runs"]["tumour_seq_type"],
               genome_build=CFG["runs"]["tumour_genome_build"],
               tumour_id=CFG["runs"]["tumour_sample_id"],
               normal_id=CFG["runs"]["normal_sample_id"],
               pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
