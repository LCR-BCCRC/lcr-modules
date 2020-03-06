#!/usr/bin/env snakemake

##### MODULES #####

from os.path  import join
from snakemake.utils import validate, min_version
from modutils import (setup_module,
                      setup_subdirs,
                      cleanup_module,
                      generate_runs,
                      symlink,
                      collapse,
                      locate_genome_bams)


##### SETUP #####

min_version("5.0.0")

CFG = setup_module(config, "manta", "1.0")
validate(CFG["samples"], schema="config/schema.yaml")
CFG = setup_subdirs(CFG, "manta", "bedpe")
CFG["runs"] = generate_runs(CFG["samples"], subgroups=[])

localrules: manta_input, manta_configure, manta_output, manta_all


##### RULES #####

rule manta_input:
    input:
        CFG["inputs"].get("sample_bam") or unpack(locate_genome_bams)
    output:
        sample_bam = join(CFG["dirs"]["input"], "{sample_id}.bam")
    run:
        symlink(input.sample_bam, output.sample_bam)
        symlink(input.sample_bam + ".bai", output.sample_bam + ".bai")


rule manta_configure:
    input:
        tumour_bam = join(CFG["dirs"]["input"], "{tumour_id}.bam"),
        normal_bam = join(CFG["dirs"]["input"], "{normal_id}.bam")
    output:
        runwf = join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}/runWorkflow.py")
    log:
        join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}/log.config.txt")
    params:
        opts   = CFG["options"]["configure"],
        fasta  = config["reference"]["genome_fasta"],
        outdir = join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        collapse("""
        configManta.py {params.opts} --referenceFasta {params.fasta} --runDir {params.outdir}
        --tumourBam {input.tumour_bam} --normalBam {input.normal_bam} > {log} 2>&1
        """)


rule manta_run:
    input:
        runwf = rules.manta_configure.output.runwf
    output:
        vcf = join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "results",
                   "variants", "somaticSV.vcf.gz")
    log:
        join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "log.manta.txt")
    params:
        opts   = CFG["options"]["manta"],
        outdir = join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}")
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
        rm -rf {params.outdir}/workspace/
        """)


rule manta_fix_vcf_ids:
    input:
        vcf  = rules.manta_run.output.vcf
    output:
        vcf = pipe(join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "results",
                        "variants", "somaticSV.with_ids.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "results", "variants", 
             "log.fix_vcf_ids.txt")
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
        vcf = pipe(join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "results",
                        "variants", "somaticSV.with_ids.with_vaf.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{tumour_id}--{normal_id}", "results", "variants", 
             "log.calc_vaf.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        collapse("""
        {input.cvaf} {input.vcf} > {output.vcf} 2> {log}
        """)


rule manta_vcf_to_bedpe:
    input:
        vcf  = rules.manta_calc_vaf.output.vcf
    output:
        bedpe = join(CFG["dirs"]["bedpe"], "{tumour_id}--{normal_id}", "somaticSV.bedpe")
    log:
        join(CFG["dirs"]["bedpe"], "{tumour_id}--{normal_id}", "log.vcf_to_bedpe.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    threads:
        CFG["threads"].get("vcf_to_bedpe") or 1
    resources: 
        mem_mb = CFG["memory"].get("vcf_to_bedpe") or 1000
    shell:
        collapse("""
        svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log}
        """)


rule manta_output:
    input:
        bedpe = rules.manta_vcf_to_bedpe.output.bedpe
    output:
        bedpe = join(CFG["dirs"]["output"], "{tumour_id}--{normal_id}.bedpe")
    run:
        symlink(input.bedpe, output.bedpe)


rule manta_all:
    input:
        vcfs = expand(rules.manta_vcf_to_bedpe.output.bedpe, zip,
                      tumour_id=CFG["runs"]["tumour_sample_id"],
                      normal_id=CFG["runs"]["normal_sample_id"])


##### CLEANUP #####

cleanup_module(CFG)

del CFG
