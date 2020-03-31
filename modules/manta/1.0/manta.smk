#!/usr/bin/env snakemake

##### MODULES #####

from os.path import join
import modutils as md


##### SETUP #####

CFG = md.setup_module(
    config = config, 
    name = "manta", 
    version = "1.0",
    subdirs = ["inputs", "manta", "bedpe", "outputs"]
)

localrules: _manta_input, _manta_configure, _manta_output, _manta_all


##### RULES #####

rule _manta_input:
    input:
        CFG["inputs"].get("sample_bam") or unpack(md.locate_bam(CFG.get("bam_directory")))
    output:
        sample_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{sample_id}.bam")
    run:
        md.symlink(input.sample_bam, output.sample_bam)
        md.symlink(input.sample_bam + ".bai", output.sample_bam + ".bai")


rule _manta_configure:
    input:
        tumour_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{tumour_id}.bam"),
        normal_bam = join(CFG["dirs"]["inputs"], "{seq_type}", "{normal_id}.bam")
    output:
        runwf = join(CFG["dirs"]["manta"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py")
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", 
             "{tumour_id}--{normal_id}--{pair_status}/log.config.txt")
    params:
        opts   = md.make_seqtype_specific(CFG["options"]["configure"]),
        fasta  = config["reference"]["genome_fasta"]
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        md.collapse("""
        configManta.py {params.opts} --referenceFasta {params.fasta} 
        --runDir "$(dirname {output.runwf})" --tumourBam {input.tumour_bam} > {log} 2>&1
        """)


rule _manta_run:
    input:
        runwf = rules._manta_configure.output.runwf
    output:
        vcf = join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{pair_status}", 
                   "results", "variants", "somaticSV.vcf.gz")
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{pair_status}", 
             "log.manta.txt")
    params:
        opts   = CFG["options"]["manta"]
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    threads:
        CFG["threads"].get("manta") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("manta") or 1000
    shell:
        md.collapse("""
        {input.runwf} {params.opts} --jobs {threads} > {log} 2>&1
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


rule _manta_fix_vcf_ids:
    input:
        vcf  = rules._manta_run.output.vcf
    output:
        vcf = pipe(join(CFG["dirs"]["manta"], "{seq_type}", 
                        "{tumour_id}--{normal_id}--{pair_status}", 
                        "results", "variants", "somaticSV.with_ids.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{pair_status}", 
             "results", "variants", "log.fix_vcf_ids.txt")
    shell:
        md.collapse("""
        gzip -dc {input.vcf}
            |
        awk 'BEGIN {{FS=OFS="\\t"}}
        $1 == "#CHROM" {{$10="{wildcards.normal_id}"; $11="{wildcards.tumour_id}"}}
        {{print $0}}' > {output.vcf} 2> {log}
        """)


rule _manta_calc_vaf:
    input:
        vcf  = rules._manta_fix_vcf_ids.output.vcf,
        cvaf = CFG["inputs"]["calc_manta_vaf"]
    output:
        vcf = pipe(join(CFG["dirs"]["manta"], "{seq_type}", 
                        "{tumour_id}--{normal_id}--{pair_status}", 
                        "results", "variants", "somaticSV.with_ids.with_vaf.vcf"))
    log:
        join(CFG["dirs"]["manta"], "{seq_type}", "{tumour_id}--{normal_id}--{pair_status}", 
             "results", "variants", "log.calc_vaf.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    shell:
        "{input.cvaf} {input.vcf} > {output.vcf} 2> {log}"


rule _manta_vcf_to_bedpe:
    input:
        vcf  = rules._manta_calc_vaf.output.vcf
    output:
        bedpe = join(CFG["dirs"]["bedpe"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{pair_status}", "somaticSV.bedpe")
    log:
        join(CFG["dirs"]["bedpe"], "{seq_type}", "{tumour_id}--{normal_id}--{pair_status}", 
             "log.vcf_to_bedpe.txt")
    conda:
        CFG["conda_envs"]["manta"] or "envs/manta.yaml"
    threads:
        CFG["threads"].get("vcf_to_bedpe") or 1
    resources: 
        mem_mb = CFG["mem_mb"].get("vcf_to_bedpe") or 1000
    shell:
        "svtools vcftobedpe -i {input.vcf} > {output.bedpe} 2> {log}"


rule _manta_output:
    input:
        bedpe = rules._manta_vcf_to_bedpe.output.bedpe
    output:
        bedpe = join(CFG["dirs"]["outputs"], "{seq_type}", 
                     "{tumour_id}--{normal_id}--{pair_status}.bedpe")
    run:
        md.symlink(input.bedpe, output.bedpe)


rule _manta_all:
    input:
        vcfs = expand(rules._manta_vcf_to_bedpe.output.bedpe, zip,
                      seq_type=CFG["paired_runs"]["tumour_seq_type"],
                      tumour_id=CFG["paired_runs"]["tumour_sample_id"],
                      normal_id=CFG["paired_runs"]["normal_sample_id"],
                      pair_status=CFG["paired_runs"]["pair_status"])


##### CLEANUP #####

md.cleanup_module(CFG)

del CFG
