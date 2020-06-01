#!/usr/bin/env snakemake


##### SETUP #####

import modutils as md

assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

CFG = md.setup_module(
    config = config, 
    name = "ginkgo", 
    version = "1.0",
    subdirs = ["inputs", "bed", "ginkgo", "outputs"],
)

localrules: 
    _ginkgo_input_bam,
    _ginkgo_link_bed_to_bins,
    _ginkgo_create_bed_list,
    _ginkgo_output,
    _ginkgo_all

##### RULES #####
# in single cell modules, sample_id and lib_id refers to patient_id and sample_id respectively on the datasheet
# for ginkgo purposes, each library = a cell

rule _ginkgo_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "{genome_build}/{sample_id}/{lib_id}.bam"
    run:
        md.symlink(input.bam, output.bam)


rule _ginkgo_bam2bed:
    input:
        bam = CFG["dirs"]["inputs"] + "{genome_build}/{sample_id}/{lib_id}.bam"
    output:
        bed = CFG["dirs"]["bed"] + "{genome_build}/{sample_id}/{lib_id}.bed.gz"
    conda:
        CFG["conda_envs"].get("bedtools") or "envs/bedtools-2.25.0.yaml"
    log:
        CFG["logs"]["bed"] + "{genome_build}/{sample_id}/{lib_id}.bedtools.stderr.log"
    threads:
        CFG["threads"].get("bam2bed") or 2
    resources:
        mem_mb = CFG["mem_mb"].get("bam2bed") or 5000
    shell:
        "bedtools bamtobed -i {input} > >(gzip > {output}) 2> {log}"


rule _ginkgo_link_bed_to_bins:
    input:
        bed = rules._ginkgo_bam2bed.output.bed
    output:
        bed = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/{lib_id}.bed.gz"
    run:
        md.symlink(input.bed, output.bed)


def create_map_dict(df = CFG["samples"]):
    sample_lib_dict = df.groupby('patient_id')['sample_id'].apply(list).to_dict()
    return sample_lib_dict
    

def get_libraries(mconfig = CFG):
    path = mconfig["dirs"]["ginkgo"]
    def _get_custom_lib(wildcards):
        d = create_map_dict()
        libs = d.get(wildcards.sample_id)
        if libs:
            return expand("{DIR}{{genome_build}}_bin{{bin}}/{{sample_id}}/{lib}.bed.gz", DIR = path, lib = libs)
        else:
            raise ValueError (f"Invalid value for patient_id : {wildcards.sample_id}")
    return _get_custom_lib

def get_lib_str(wildcards):
    d = create_map_dict()
    libs = d.get(wildcards.sample_id)
    print(len(libs))
    return "\|".join(libs)


rule _ginkgo_create_bed_list:
    input:
        bed = get_libraries(CFG),
        metrics = CFG["inputs"]["metrics"]
    output:
        bed_list = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/bed_files.txt"
    params:
        lib_str = get_lib_str,
        bedDir = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}",
        opts = CFG["options"]["bed_threshold"]
    shell:
        md.as_one_line("""
        grep \"{params.lib_str}\" {input.metrics} | awk -v P=`pwd` 'NR > 1 {{if ($2 >= {params.opts}) {{ print P \"/{params.bedDir}/\" $1 \".bed.gz\" }}}}' > {output.bed_list}
        """)


rule _ginkgo_run:
    input:
        bed_list = rules._ginkgo_create_bed_list.output.bed_list 
        # get_libraries(CFG) #expand("{DIR}{{genome_build}}_bin{{bin}}/{{sample_id}}/{lib}.bed.gz", DIR = CFG["dirs"]["ginkgo"], lib = get_libraries)
    output:
        sc = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}/SegCopy",
        stamp = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.done"
    conda:
        CFG["conda_envs"].get("ginkgo") or "envs/ginkgo_dep.yaml"
    log: 
        stdout = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}.stderr.log"
    params:
        ginkgo = CFG["inputs"]["script"],
        bedDir = CFG["dirs"]["ginkgo"] + "{genome_build}_bin{bin}/{sample_id}",
        genome = "{genome_build}",
        binning = CFG["options"]["binMeth"],
        clustDist = CFG["options"]["distMeth"],
        clustLink = CFG["options"]["clustMeth"],
        opts = CFG["options"]["flags"]
    threads: 
        CFG["threads"].get("ginkgo") or 4
    resources:
        mem_mb = CFG["mem_mb"].get("ginkgo") or 8000
    shell:
        md.as_one_line("""
        {params.ginkgo} 
        --input {params.bedDir} 
        --genome {params.genome} 
        --binning {params.binning}
        --cells {input.bed_list}
        {params.clustDist} {params.clustLink}
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        && touch {output.stamp}
        """)


rule _ginkgo_output:
    input:
        sc = rules._ginkgo_run.output.sc
    output:
        sc = CFG["dirs"]["outputs"] + "{genome_build}_bin{bin}/{sample_id}_SegCopy"
    run:
        md.symlink(input.sc, output.sc)

sub_df = CFG["samples"][["patient_id", "genome_build"]].drop_duplicates()

rule _ginkgo_all:
    input:
        expand(expand("{{dir}}{genome_build}_bin{{bin}}/{sample_id}_SegCopy", zip,
               genome_build = sub_df["genome_build"],
               sample_id = sub_df["patient_id"]),
            dir = CFG["dirs"]["outputs"],
            bin = CFG["inputs"]["bins"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
