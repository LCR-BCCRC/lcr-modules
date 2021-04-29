#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Nicole Thomas
# Module Author:    Jasper Wong
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["gatk_rnaseq"]`
CFG = op.setup_module(
    name = "gatk_rnaseq",
    version = "1.0",
    subdirectories = ["inputs", "gatk_splitntrim", "base_recal_report", "gatk_applybqsr", "gatk_variant_calling", "merge_vcfs", "gatk_variant_filtration", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gatk_rnaseq_input_bam,
    _gatk_rnaseq_output_vcf,
    _gatk_rnaseq_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gatk_rnaseq_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)

rule _gatk_splitntrim:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["gatk_splitntrim"] +  "bam/{seq_type}--{genome_build}--{pair_status}/{sample_id}.split_reassign_mq.bam")
    log:
        stdout = CFG["logs"]["gatk_splitntrim"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_splitntrim.stdout.log",
        stderr = CFG["logs"]["gatk_splitntrim"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_splitntrim.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb / 1000)
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads:
        CFG["threads"]["gatk_splitntrim"]
    resources:
        **CFG["resources"]["gatk_splitntrim"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}G" SplitNCigarReads -fixNDN TRUE -RF GoodCigarReadFilter
        -R {input.fasta} -I {input.bam} -O {output.bam}
        > {log.stdout} 2> {log.stderr}
        """)

rule _gatk_base_recalibration:
    input:
        bam = rules._gatk_splitntrim.output,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        table = CFG["dirs"]["base_recal_report"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.recalibration_report.grp"
    log:
        stdout = CFG["logs"]["base_recal_report"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["base_recal_report"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        dbsnp = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        mem_mb = lambda wildcards, resources: int(resources.mem_mb / 1000)
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_base_recalibration"]
    resources:
        **CFG["resources"]["gatk_base_recalibration"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}G" BaseRecalibrator 
        -R {input.fasta} -I {input.bam} -known-sites {params.dbsnp} -known-sites {params.gnomad} -O {output.table} > {log.stdout} 2> {log.stderr}
        """)
        
rule _gatk_applybqsr:
    input:
        bam = rules._gatk_splitntrim.output.bam,
        table = rules._gatk_base_recalibration.output.table,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        bam = temp(CFG["dirs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.recalibrated.bam"),
        bai = temp(CFG["dirs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.recalibrated.bai")
    log:
        stdout = CFG["logs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_applybqsr"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb / 1000)
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_applybqsr"]
    resources:
        **CFG["resources"]["gatk_applybqsr"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}G" ApplyBQSR 
        -R {input.fasta} -I {input.bam} -bqsr-recal-file {input.table} -O {output.bam} > {log.stdout} 2> {log.stderr}
        """)
        
# Symlink chromosomes used for parallelization
checkpoint _gatk_rnaseq_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.chrs, output.chrs)
        
        
rule _gatk_variant_calling:
    input:
        bam = rules._gatk_applybqsr.output.bam,
        bai = rules._gatk_applybqsr.output.bai,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = temp(CFG["dirs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.{chrom}.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.{chrom}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_variant_calling"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.{chrom}.gatk_base_recal.stderr.log"
    params:
        # java_opts = "-Xmx{0}G".format(int(CFG["resources"]["gatk_variant_calling"]["mem_mb"]) // 1000),
        mem_mb = lambda wildcards, resources: int(resources.mem_mb / 1000),
        vcf = reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz"),
        min_conf_thres = CFG["options"]["gatk_variant_calling"]["min_conf_thres"]
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_variant_calling"]
    resources:
        **CFG["resources"]["gatk_variant_calling"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}G" HaplotypeCaller 
        -R {input.fasta} -I {input.bam} -dont-use-soft-clipped-bases -stand-call-conf {params.min_conf_thres} -L {wildcards.chrom}
        -dbsnp {params.vcf} -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)
        
        
def _gatk_rnaseq_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["gatk_rnaseq"]
    chrs = checkpoints._gatk_rnaseq_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["gatk_variant_calling"] + "{{seq_type}}--{{genome_build}}--{{pair_status}}/{{sample_id}}.{chrom}.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _gatk_rnaseq_get_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["gatk_rnaseq"]
    chrs = checkpoints._gatk_rnaseq_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["gatk_variant_calling"] + "{{seq_type}}--{{genome_build}}--{{pair_status}}/{{sample_id}}.{chrom}.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


# Merge chromosome gatk_rnaseq VCFs from the same sample
rule _gatk_rnaseq_merge_vcfs:
    input:
        vcf = _gatk_rnaseq_get_chr_vcfs,
        tbi = _gatk_rnaseq_get_chr_tbis
    output:
        vcf = temp(CFG["dirs"]["merge_vcfs"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["merge_vcfs"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.output.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["merge_vcfs"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.merge_vcfs_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["merge_vcfs"]
    resources:
        **CFG["resources"]["merge_vcfs"]
    params: 
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -m {params.mem_mb}M -O z -o {output.vcf} 2>> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)
        
        
rule _gatk_variant_filtration:
    input:
        vcf = rules._gatk_rnaseq_merge_vcfs.output.vcf,
        tbi = rules._gatk_rnaseq_merge_vcfs.output.tbi,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.filtered.vcf.gz",
        tbi = CFG["dirs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.filtered.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stdout.log",
        stderr = CFG["logs"]["gatk_variant_filtration"] + "{seq_type}--{genome_build}--{pair_status}/{sample_id}.gatk_base_recal.stderr.log"
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb / 1000)
    conda:
        CFG["conda_envs"]["gatk_rnaseq"]
    threads: CFG["threads"]["gatk_variant_filtration"]
    resources:
        **CFG["resources"]["gatk_variant_filtration"]
    shell:
        # FS - phred score with strand bias > 30; QD - Variant conf/qual by depth > 5; DP - read depth > 5
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}G" VariantFiltration -R {input.fasta} -V {input.vcf} 
        -window 35 -cluster-size 3  
        -filter-expression "FS > 30.0" -filter-name FS
        -filter-expression "QD > 5.0" -filter-name QD 
        -filter-expression "DP > 5.0" -filter-name DP 
        -O {output.vcf} > {log.stdout} 2> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gatk_rnaseq_output_vcf:
    input:
        vcf = str(rules._gatk_variant_filtration.output.vcf),
        tbi = str(rules._gatk_variant_filtration.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}--{pair_status}/{sample_id}.output.filt.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}--{pair_status}/{sample_id}.output.filt.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


# Generates the target sentinels for each run, which generate the symlinks
rule _gatk_rnaseq_all:
    input:
        expand(
            [
                str(rules._gatk_rnaseq_output_vcf.output.vcf),
                str(rules._gatk_rnaseq_output_vcf.output.tbi)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            sample_id=CFG["runs"]["tumour_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
