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
    subdirectories = ["inputs", "lofreq", "combined", "filtered", "outputs"],
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
        vcf_snvs_filtered = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.snvs.vcf.gz",
        vcf_indels_filtered = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.indels.vcf.gz",
        vcf_snvs_all = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.snvs.vcf.gz",
        vcf_indels_all = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.indels.vcf.gz"
    log:
        stdout = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stdout.log",
        stderr = CFG["logs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq.stderr.log"
    params:
        opts = CFG["options"]["lofreq"],
        regions = op.switch_on_wildcard("seq_type", CFG["switches"]["regions_bed"]),
        rm_files = CFG["dirs"]["lofreq"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/normal_relaxed.vcf.gz"
    conda:
        CFG["conda_envs"]["lofreq"]
    threads:
        CFG["threads"]["lofreq"]
    resources:
        **CFG["resources"]["lofreq"]
    shell:
        op.as_one_line("""
        lofreq somatic {params.opts} --threads {threads} -t {input.tumour_bam} -n {input.normal_bam}
        -f {input.fasta} -o $(dirname {output.vcf_snvs_filtered})/ -d {input.dbsnp} {params.regions} 
        > {log.stdout} 2> {log.stderr}
            &&
        rm {params.rm_files}
        """)


rule _lofreq_combine_vcf:
    input:
        vcf_all = expand(CFG["dirs"]["lofreq"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/somatic_final.{var_type}.vcf.gz",
                    var_type = ["indels", "snvs"]),
        vcf_all_filtered = expand(CFG["dirs"]["lofreq"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/somatic_final_minus-dbsnp.{var_type}.vcf.gz",
                    var_type = ["indels", "snvs"])
    output:
        vcf_all = temp(CFG["dirs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.vcf.gz"),
        vcf_all_filtered = temp(CFG["dirs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.vcf.gz"),
    resources:
        **CFG["resources"]["bcftools_sort"]
    conda:
        CFG["conda_envs"]["bcftools"]
    log:
        stdout_all = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final.combined.stdout.log",
        stderr_all = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final.combined.stderr.log",
        stdout_all_filtered = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final_minus-dbsnp.combined.stdout.log",
        stderr_all_filtered = CFG["logs"]["combined"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lofreq_final_minus-dbsnp.combined.stderr.log"
    shell:
                op.as_one_line("""
        bcftools concat -a {input.vcf_all} |
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf_all}
        > {log.stdout_all} 2> {log.stderr_all}
            &&
        bcftools concat -a {input.vcf_all_filtered} |
        bcftools sort --max-mem {resources.mem_mb}M -Oz -o {output.vcf_all_filtered}
        > {log.stdout_all_filtered} 2> {log.stderr_all_filtered}
        """)


rule _lofreq_filter_vcf:
    input:
        vcf_all = rules._lofreq_combine_vcf.output.vcf_all,
        vcf_all_filtered = rules._lofreq_combine_vcf.output.vcf_all_filtered,
        lofreq_filter = CFG["inputs"]["lofreq_filter"]
    output:
        vcf_all_clean = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.filtered.vcf.gz",
        vcf_all_filtered_clean = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.filtered.vcf.gz",
        vcf_all_clean_tbi = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final.combined.filtered.vcf.gz.tbi",
        vcf_all_filtered_clean_tbi = CFG["dirs"]["filtered"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somatic_final_minus-dbsnp.combined.filtered.vcf.gz.tbi"
    resources:
        **CFG["resources"]["bcftools_sort"]
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bash {input.lofreq_filter} {input.vcf_all} | bgzip > {output.vcf_all_clean}
          && tabix -p vcf {output.vcf_all_clean}
              &&
        bash {input.lofreq_filter} {input.vcf_all_filtered} | bgzip > {output.vcf_all_filtered_clean}
          && tabix -p vcf {output.vcf_all_filtered_clean}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lofreq_output_vcf:
    input:
        vcf_all = rules._lofreq_filter_vcf.output.vcf_all_clean,
        vcf_all_filtered = rules._lofreq_filter_vcf.output.vcf_all_filtered_clean
    output:
        vcf_all = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.snvs.vcf.gz",
        vcf_all_filtered = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_minus-dbsnp.{tool}.snvs.vcf.gz"
    run:
        op.relative_symlink(input.vcf_all, output.vcf_all)
        op.relative_symlink(input.vcf_all + ".tbi", output.vcf_all + ".tbi")
        op.relative_symlink(input.vcf_all_filtered, output.vcf_all_filtered)
        op.relative_symlink(input.vcf_all_filtered + ".tbi", output.vcf_all_filtered + ".tbi")


# Generates the target sentinels for each run, which generate the symlinks
rule _lofreq_all:
    input:
        expand(
            [
                rules._lofreq_output_vcf.output.vcf_all,
                rules._lofreq_output_vcf.output.vcf_all_filtered,
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            tool="lofreq")


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
