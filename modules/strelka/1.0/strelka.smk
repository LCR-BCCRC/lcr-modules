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
    subdirectories = ["inputs", "strelka", "filtered", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _strelka_input_bam,
    _strelka_configure,
    _strelka_unzip,
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
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
rule _strelka_configure:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        runwf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/runWorkflow.py"
    log:
        stdout = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stdout.log",
        stderr = CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/strelka_configure.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["configure"]),
        indels = "--indelCandidates=" + CFG["inputs"].get("candidate_indel")
    conda:
        CFG["conda_envs"]["strelka"]
    shell:
        op.as_one_line("""
        configureStrelkaSomaticWorkflow.py 
        --normalBam={input.normal_bam}
        --tumorBam={input.tumour_bam}
        --referenceFasta={input.fasta}
        --runDir=$(dirname {params.outDIR})
        {params.indels}
        {params.opts} 
        > {log.stdout} 2> {log.stderr}
        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _strelka_run:
    input:
        runwf = rules._strelka_configure.output.runwf
    output:
        snp = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.snvs.vcf.gz",
        indels = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.indels.vcf.gz"
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
        mem_mb = CFG["mem_mb"]["strelka"]
    shell:
        op.as_one_line("""
        {input.runwf} -j {threads} 
        {params.opts} 
        >{log.stdout} 2> {log.stderr}
        """)


rule _strelka_unzip:
    input:
        vcf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.{var_type}.vcf.gz"
    output:
        vcf = CFG["dirs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/somatic.{var_type}.vcf"
    log:
        CFG["logs"]["strelka"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}_unzip.stderr.log"
    shell:
        "zcat {input.vcf} > {output.vcf} 2> {log}"


rule _strelka_filter:
    input:
        vcf = rules._strelka_unzip.output.vcf
    output:
        vcf = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/results/variants/{var_type}.filt.vcf"
    shell:
        op.as_one_line("""
        grep '^#' {input.vcf} > {output.vcf} && 
        awk '{{ if ($7 == \"PASS\") {{ print $0 }} }}' {input.vcf} >> {output.vcf}
        """)



# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _strelka_output_vcf:
    input:
        vcf = rules._strelka_filter.output.vcf
    output:
        vcf = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{var_type}.filt.vcf"
    run:
        op.relative_symlink(input, output)


# Generates the target sentinels for each run, which generate the symlinks
rule _strelka_all:
    input:
        expand(expand("{dir}{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{{var_type}}.filt.vcf", zip,
                dir = CFG["dirs"]["outputs"],
                seq_type = CFG["runs"]["tumour_seq_type"],
                genome_build = CFG["runs"]["tumour_genome_build"],
                tumour_id = CFG["runs"]["tumour_sample_id"],
                normal_id = CFG["runs"]["normal_sample_id"],
                pair_status = CFG["runs"]["pair_status"]),
            var_type = ['indels', 'snvs'])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
