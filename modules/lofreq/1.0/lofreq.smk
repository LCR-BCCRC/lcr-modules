#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Bruno Grande
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["lofreq"]`
CFG = op.setup_module(
    name = "lofreq",
    version = "1.0",
    subdirectories = ["inputs", "lofreq", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _lofreq_input_bam,
    _lofreq_output_vcf,
    _lofreq_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _lofreq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


# Run LoFreq in somatic variant calling mode
rule _lofreq_run:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")
    output:
        vcf_snvs = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.snvs.vcf.gz",
        vcf_indels = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.indels.vcf.gz"
    log:
        stdout = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stdout.log",
        stderr = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stderr.log"
    params:
        opts = CFG["options"]["lofreq"],
        regions = op.switch_on_wildcard("seq_type", CFG["switches"]["regions_bed"])
    conda:
        CFG["conda_envs"]["lofreq"]
    threads:
        CFG["threads"]["lofreq"]
    resources:
        mem_mb = CFG["mem_mb"]["lofreq"],
        bam =1
    shell:
        op.as_one_line("""
        lofreq somatic {params.opts} --threads {threads} -t {input.tumour_bam} -n {input.normal_bam}
        -f {input.fasta} -o $(dirname {output.vcf_snvs})/ -d {input.dbsnp} {params.regions} 
        > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lofreq_output_vcf:
    input:
        vcf_snvs = rules._lofreq_run.output.vcf_snvs,
        vcf_indels = rules._lofreq_run.output.vcf_indels
    output:
        vcf_snvs = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.snvs.vcf.gz",
        vcf_indels = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.indels.vcf.gz"
    run:
        op.relative_symlink(input.vcf_snvs, output.vcf_snvs)
        op.relative_symlink(input.vcf_snvs + ".tbi", output.vcf_snvs + ".tbi")
        op.relative_symlink(input.vcf_indels, output.vcf_indels)
        op.relative_symlink(input.vcf_indels + ".tbi", output.vcf_indels + ".tbi")


# Generates the target sentinels for each run, which generate the symlinks
rule _lofreq_all:
    input:
        expand(
            [
                rules._lofreq_output_vcf.output.vcf_snvs,
                rules._lofreq_output_vcf.output.vcf_indels,
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
