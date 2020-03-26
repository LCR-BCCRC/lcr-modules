#!/usr/bin/env snakemake

##### MODULES #####

from os.path  import join
from modutils import (setup_module,
                      cleanup_module,
                      symlink,
                      collapse,
                      locate_genome_bams)


##### SETUP #####

CFG = setup_module(
    config = config, 
    name = "manta", 
    version = "1.0",
    subdirs = ["manta", "bedpe"],
    return_unpaired=True
)

localrules: manta_input, manta_configure, manta_output, manta_all


##### RULES #####

rule manta_input:
    input:
        CFG["inputs"].get("sample_bam") or unpack(locate_genome_bams)
    output:
        sample_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{sample_id}.bam")
    run:
        symlink(input.sample_bam, output.sample_bam)
        symlink(input.sample_bam + ".bai", output.sample_bam + ".bai")


rule manta_configure:
    input:
        tumour_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{tumour_id}.bam"),
        normal_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{normal_id}.bam")
    output:
        runwf = join(CFG["dirs"]["manta"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{is_matched}/runWorkflow.py")
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", 
             "{tumour_id}--{normal_id}--{is_matched}/log.config.txt")
    params:
        opts   = CFG["options"]["configure"],
        fasta  = config["reference"]["genome_fasta"]
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        collapse("""
        configManta.py {params.opts} --referenceFasta {params.fasta} 
        --runDir "$(dirname {output.runwf})" --tumourBam {input.tumour_bam} > {log} 2>&1
        """)


rule manta_run:
    input:
        runwf = rules.manta_configure.output.runwf
    output:
        vcf = join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", 
                   "results", "variants", "somaticSV.vcf.gz")
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", 
             "log.manta.txt")
    params:
        opts   = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    threads:
        CFG["threads"].get("manta") or 1
    resources: 
        mem_mb = CFG["memory"].get("manta") or 1000
    shell:
        collapse("""
        {input.runwf} {params.opts} --jobs {threads} > {log} 2>&1
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


rule manta_fix_vcf_ids:
    input:
        vcf  = rules.manta_run.output.vcf
    output:
        vcf = pipe(join(CFG["dirs"]["manta"], "{seq_type}", 
                        "{tumour_id}--{normal_id}--{is_matched}", 
                        "results", "variants", "somaticSV.with_ids.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", 
             "results", "variants", "log.fix_vcf_ids.txt")
    shell:
        collapse("""
        gzip -dc {input.vcf}
            |
        awk 'BEGIN {{FS=OFS="\\t"}}
        $1 == "#CHROM" {{$10="{wildcards.normal_id}"; $11="{wildcards.tumour_id}"}}
        {{print $0}}' > {output.vcf} 2> {log}
        """)


rule manta_calc_vaf:
    input:
        vcf  = rules.manta_fix_vcf_ids.output.vcf,
        cvaf = CFG["inputs"]["calc_manta_vaf"]
    output:
        vcf = pipe(join(CFG["dirs"]["manta"], "{seq_type}", 
                        "{tumour_id}--{normal_id}--{is_matched}", 
                        "results", "variants", "somaticSV.with_ids.with_vaf.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", 
             "results", "variants", "log.calc_vaf.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        "{input.cvaf} {input.vcf} > {output.vcf} 2> {log}"


rule manta_vcf_to_bedpe:
    input:
        vcf  = rules.manta_calc_vaf.output.vcf
    output:
        bedpe = join(CFG["dirs"]["bedpe"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{is_matched}", "somaticSV.bedpe")
    log:
        join(CFG["dirs"]["bedpe"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", 
             "log.vcf_to_bedpe.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    threads:
        CFG["threads"].get("vcf_to_bedpe") or 1
    resources: 
        mem_mb = CFG["memory"].get("vcf_to_bedpe") or 1000
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log}"


rule manta_output:
    input:
        bedpe = rules.manta_vcf_to_bedpe.output.bedpe
    output:
        bedpe = join(CFG["dirs"]["outputs"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{is_matched}.bedpe")
    run:
        symlink(input.bedpe, output.bedpe)


rule manta_all:
    input:
        vcfs = expand(rules.manta_vcf_to_bedpe.output.bedpe, zip,
                      seq_type=CFG["runs"]["tumour_seq_type"],
                      tumour_id=CFG["runs"]["tumour_sample_id"],
                      normal_id=CFG["runs"]["normal_sample_id"],
                      is_matched=CFG["runs"]["is_matched"])


##### CLEANUP #####

cleanup_module(CFG)

del CFG
