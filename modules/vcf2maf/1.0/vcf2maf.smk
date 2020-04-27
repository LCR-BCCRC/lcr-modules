#!/usr/bin/env snakemake

##### MODULES #####

from os.path  import join
from modutils import (setup_module,
                      cleanup_module,
                      symlink,
                      collapse)


##### SETUP #####

CFG = setup_module(
    config = config, 
    name = "vcf2maf", 
    version = "1.0",
    subdirs = ["vcf2maf"],
    return_unpaired=True
)

localrules: vcf2maf_input, vcf2maf_output, vcf2maf_all


##### RULES #####

rule vcf2maf_input:
    input:
        vcf = CFG["inputs"].get("vcf")
    output:
        vcf = join(CFG["dirs"]["inputs"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.vcf")
    run:
        symlink(input.vcf, output.vcf)


rule vcf2maf_run:
    input:
        vcf = rules.vcf2maf_input.output.vcf
    output:
        maf = join(CFG["dirs"]["vcf2maf"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.maf")
    log:
        join(CFG["dirs"]["vcf2maf"], "{seq_type}", 
             "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.log.txt")
    params:
        opts = CFG["options"]["vcf2maf"],
        fasta  = config["reference"]["genome_fasta"],
        vep = CFG["inputs"]["vep"]
    conda:
        CFG["conda_envs"]["vcf2maf"] or "envs/vcf2maf.yaml"
    threads:
        CFG["threads"].get("vcf2maf") or 1
    resources: 
        mem_mb = CFG["memory"].get("vcf2maf") or 4000
    shell:
        """
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --ref-fasta {params.fasta} --vep-data {params.vep} --vep-path $vepPATH {params.opts} 2> {log}
        """
        

rule vcf2maf_output:
    input:
        maf = rules.vcf2maf_run.output.maf
    output:
        maf = join(CFG["dirs"]["outputs"], "{seq_type}", "{tumour_id}--{normal_id}--{is_matched}", "{caller}", "{var_type}.maf"),
    run:
        symlink(input.maf, output.maf)

rule vcf2maf_all:
    input:
        vcfs = expand(expand("{dir}/{seq_type}/{tumour_id}--{normal_id}--{is_matched}/{{caller}}/{{var_type}}.maf", zip,
                    dir = CFG["dirs"]["outputs"],
                    seq_type = CFG["runs"]["tumour_seq_type"],
                    tumour_id = CFG["runs"]["tumour_sample_id"],
                    normal_id = CFG["runs"]["normal_sample_id"],
                    is_matched = CFG["runs"]["is_matched"]),
                var_type = ['indels', 'snvs'], caller = VCFcaller)


##### CLEANUP #####

cleanup_module(CFG)

del CFG
