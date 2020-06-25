#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["strelka"]`
CFG = op.setup_module(
    name = "strelka",
    version = "1.0",
    subdirectories = ["inputs", "chrom_bed", "strelka", "filtered", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _strelka_input_bam,
    _strelka_input_vcf,
    _strelka_dummy_vcf,
    _strelka_index_bed,
    _strelka_configure_paired,
    _strelka_configure_unpaired,
    _strelka_decompress_vcf,
    _strelka_filter,
    _strelka_output_vcf,
    _strelka_all,


##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _strelka_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


rule _strelka_input_vcf:
    input:
        vcf = CFG["inputs"]["candidate_small_indels_vcf"],
        tbi = CFG["inputs"]["candidate_small_indels_tbi"]
    output:
        vcf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.candidateSmallIndels.vcf.gz",
        tbi = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.candidateSmallIndels.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


rule _strelka_dummy_vcf:
    output:
        touch(CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/vcf/{tumour_id}--{normal_id}--{pair_status}.dummy")


# bgzip-compress and tabix-index the BED file to meet strelka requirement
rule _strelka_index_bed:
    input:
        bed = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.bed")
    output:
        bedz = CFG["dirs"]["chrom_bed"] + "{genome_build}.main_chroms.bed.gz"
    conda:
        CFG["conda_envs"]["tabix"]
    shell:
        op.as_one_line("""
        bgzip -c {input.bed} > {output.bedz}
            &&
        tabix {output.bedz}
        """)


def _get_indel_cli_arg(wildcards, input):
    if "dummy" in input.indels:
        param = ""
    else: 
        param = f"--indelCandidates={input.indels}"
    return param


rule _strelka_configure_paired: # Somatic
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bedz = rules._strelka_index_bed.output.bedz,
        indels = rules._strelka_input_vcf.output.vcf if CFG["inputs"]["candidate_small_indels_vcf"] else rules._strelka_dummy_vcf.output
    output:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stderr.log"
    params:
        indel_arg = _get_indel_cli_arg,
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
    wildcard_constraints:
        pair_status = "matched|unmatched"
    conda:
        CFG["conda_envs"]["strelka"]
    shell:
        op.as_one_line("""
        configureStrelkaSomaticWorkflow.py 
        --normalBam={input.normal_bam}
        --tumorBam={input.tumour_bam}
        --referenceFasta={input.fasta}
        --callRegions={input.bedz}
        --runDir=$(dirname {output.runwf})
        {params.indel_arg}
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _strelka_configure_unpaired: # germline
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        bedz = rules._strelka_index_bed.output.bedz
    output:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"])
    message:
        "WARNING: {wildcards.seq_type} sample ({wildcards.sample}) is processed using Strelka Germline workflow. Ensure pairing config for capture is set to run_unpaired_tumours_with: 'unmatched_normal' to run Strelka Somatic workflow"
    wildcard_constraints:
        pair_status = "no_normal"
    conda:
        CFG["conda_envs"]["strelka"]
    shell:
        op.as_one_line("""
        configureStrelkaGermlineWorkflow.py 
        --bam={input.tumour_bam}
        --referenceFasta={input.fasta}
        --callRegions={input.bedz}
        --runDir=$(dirname {output.runwf})
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)


rule _strelka_run:
    input:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    output:
        vcf_dir = directory(CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/")
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_run.stderr.log"
    params:
        opts = CFG["options"]["strelka"]
    conda:
        CFG["conda_envs"]["strelka"]
    threads:
        CFG["threads"]["strelka"]
    resources: 
        mem_mb = op.retry(CFG["mem_mb"]["strelka"], 2)
    shell:
        op.as_one_line("""
        {input.runwf} {params.opts} --jobs {threads} > {log.stdout} 2> {log.stderr}
            &&
        rm -rf "$(dirname {input.runwf})/workspace/"
        """)


rule _strelka_decompress_vcf:
    input:
        vcf_dir = rules._strelka_run.output.vcf_dir
    output:
        vcf = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}.vcf"
    shell:
        "gzip -dc {input.vcf_dir}{wildcards.var_type}.vcf.gz > {output.vcf}"


rule _strelka_filter:
    input:
        vcf = rules._strelka_decompress_vcf.output.vcf
    output:
        vcf = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}.passed.vcf"
    shell:
        op.as_one_line("""
        awk 'BEGIN {{FS=OFS="\t"}} $0 ~ /^#/ || $7 == "PASS"' {input.vcf} > {output.vcf}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _strelka_output_vcf:
    input:
        vcf = rules._strelka_decompress_vcf.output.vcf,
        vcf_p = rules._strelka_filter.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}.vcf",
        vcf_p = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}.passed.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.vcf_p, output.vcf_p)


def _get_strelka_output(wildcards):
    CFG = config["lcr-modules"]["strelka"]

    if wildcards.pair_status == "no_normal":
        vcf_files = "variants"
    else:
        vcf_files = ["somatic.snvs", "somatic.indels"]
    vcf = expand(rules._strelka_output_vcf.output.vcf, var_type = vcf_files, **wildcards)
    return vcf


rule _strelka_dispatch:
    input: 
        vcf = _get_strelka_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _strelka_all:
    input:
        expand(rules._strelka_dispatch.output.dispatched, zip,
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"])

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
