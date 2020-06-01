#!/usr/bin/env snakemake


##### SETUP #####
import modutils as md

# Make sure the `CFG` variable doesn"t exist yet
assert "CFG" not in locals(), "`CFG` is a reserved variable for lcr-modules."

# Setup module and store module-specific configuration in `CFG`.
# `CFG` corresponds to `config["lcr-modules"]["cellranger"]`

localrules: 
    _cellranger_input,
    _cellranger_output,
    _cellranger_all


##### RULES #####
rule _cellranger_create_samplesheet:
    output:
        ss = CFG["dirs"]["samplesheet"] + "{chip_id}_samplesheet.csv"
    run:
        df = CFG["samples"]
        ss = df[["lane", "sample_id", "index"]]
        ss.columns = ["lane", "sample", "index"]
        ss.to_csv({output.ss}, sep = ",", index = False)

def get_completion_files(raw_dir = CFG["dirs"]["raw"], suffix = ["RTAComplete*", "RunInfo*", "RunParameters*"]):
    def get_custom_files(wildcards):
        path = raw_dir + "/*" + pattern + "*"
        runs = glob.glob(path)
        file = []
        for f in suffix:
            file += glob.glob(runs + "/" + f)
        return file
    return get_custom_files

rule _cellranger_mkfastq:
    input:
        check = get_completion_files
        ss = rules._cellranger_create_samplesheet
    output:
        stamp = CFG["dirs"]["mkfastq"] + "stamps/{run_id}_mkfastq.stamp"
    log:
        stdout = CFG["logs"]["mkfastq"] + "{run_id}/mkfastq.stdout.log",
        stderr = CFG["logs"]["mkfastq"] + "{run_id}/mkfastq.stderr.log"
    params:
        opts = CFG["options"]["mkfastq"],
        ss = 
    conda:
        CFG["conda_envs"].get("cellranger") or "envs/cellranger.yaml"
    threads:
        CFG["threads"]["step_1"]
    resources: 
        mem_mb = CFG["mem_mb"]["step_1"]
    shell:
        md.as_one_line("""
        <TODO> {params.opts} --bam {input.bam} --ref-fasta {params.fasta} 
        --output {output.vcf} --threads {threads} > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _cellranger_step_2:
    input:
        vcf = rules._cellranger_step_1.output.vcf
    output:
        vcf = CFG["dirs"]["cellranger"] + "{seq_type}--{genome_build}/{sample_id}/variants.filt.vcf"
    log:
        stderr = CFG["logs"]["cellranger"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "gzip -dc {input.vcf} | grep {params.opts} > {output.vcf} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under "99-outputs/").
# TODO: Update output file (and if applicable, add one rule for each output file)
rule _cellranger_output_vcf:
    input:
        vcf = rules._cellranger_step_2.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{sample_id}.variants.filt.vcf"
    run:
        md.symlink(input.vcf, output.vcf)


# Generates the target sentinels for each run, which generate the symlinks
# TODO: Update to ask for the output of every `_cellranger_output_*` rule
rule _cellranger_all:
    input:
        expand(rules._cellranger_output.output.vcf, 
               zip,  # Run expand() with zip(), not product()
               seq_type=CFG["samples"]["seq_type"],
               genome_build=CFG["samples"]["genome_build"],
               sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk (including the samples and runs)
md.cleanup_module(CFG)

# Delete any local variables to avoid interfering with other code
del CFG
