#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Lauren Chong
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["varscan"]`
CFG = op.setup_module(
    name = "varscan",
    version = "1.0",
    subdirectories = ["inputs", "mpileup", "varscan", "maf", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _varscan_input_bam,
    _varscan_input_chroms,
    _varscan_reheader_vcf,
    _varscan_combine_vcf,
    _varscan_symlink_maf,
    _varscan_output_vcf,
    _varscan_output_maf,
    _varscan_dispatch,
    _varscan_all,

ruleorder: _varscan_reheader_vcf > _varscan_combine_vcf

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _varscan_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bam + ".bai")


# Pulls in list of chromosomes for the genome builds
checkpoint _varscan_input_chroms:
    input:
        txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.txt, output.txt)

# generate mpileups for tumour and normal bams separately. 
# If we parallelize this by chromosome we will need 2 * 2 threads per chromosome but this should be a lot more efficient
rule _varscan_bam2mpu:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        mpu = temp(CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup")
    log:
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{sample_id}.bam2mpu.{chrom}.stderr.log"
    params:
        opts = CFG["options"]["mpileup"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        1    #hardcoded because samtools mpileup does not support more than one thread.
    resources:
        mem_mb = CFG["mem_mb"]["mpileup"]
    shell:
        op.as_one_line("""
        samtools mpileup {params.opts}
        -r {wildcards.chrom}
        -f {input.fasta} {input.bam}
        > {output.mpu}
        2> {log.stderr}
        """)


rule _varscan_somatic:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.{chrom}.mpileup",
        normal_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{normal_id}.{chrom}.mpileup"
    output:
        snp = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_snp.vcf"),
        indel = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_indel.vcf")
    wildcard_constraints:
        pair_status = "matched|unmatched"
    log:
        stdout = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stdout.log",
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["somatic"])
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["somatic"]   #this seems to rarely exceed 300% due to samtools and/or I/O restrictions
    resources:
        mem_mb = CFG["mem_mb"]["somatic"]
    shell:
        op.as_one_line("""
        varscan somatic 
        {input.normal_mpu} {input.tumour_mpu} 
        --output-snp {output.snp} --output-indel {output.indel}
        {params.opts}
        > {log.stdout} 2> {log.stderr} 
        || true
        """)

rule _varscan_unpaired:
    input:
        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.{chrom}.mpileup"
    output:
        vcf = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf")
    wildcard_constraints:
        pair_status = "no_normal"
    log:
        stderr = CFG["logs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_{vcf_name}.stderr.log"
    params:
        opts = op.switch_on_wildcard("seq_type", CFG["options"]["unpaired"])
    conda:
        CFG["conda_envs"]["varscan"]
    threads:
        CFG["threads"]["unpaired"]
    resources:
        mem_mb = CFG["mem_mb"]["unpaired"]
    shell:
        op.as_one_line("""
        varscan mpileup2{wildcards.vcf_name} 
        {input.tumour_mpu} 
        {params.opts}
        > {output.vcf} 2> {log.stderr}
        """)


rule _varscan_reheader_vcf:
    input:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf",
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}_{vcf_name}.vcf.gz"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        contig=$( awk '{{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}}' {input.fai})
            &&
        awk -v var="$contig" '/^#CHROM/ {{ print(var); }} {{print;}}' {input.vcf} > {output.vcf}
        """)

        #bcftools reheader -h {input.header} {input.vcf} | bcftools view -l 1 -o {output.vcf}


def _varscan_request_chrom_vcf(wildcards):
    CFG = config["lcr-modules"]["varscan"]
    with open(checkpoints._varscan_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    vcf_files = expand(str(rules._varscan_reheader_vcf.output.vcf),
        chrom = mains_chroms, **wildcards
    )
    return vcf_files


rule _varscan_combine_vcf:
    input:
        vcf = _varscan_request_chrom_vcf
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pass.somatic.{vcf_name}.vcf",
        vcf_gz = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pass.somatic.{vcf_name}.vcf.gz"
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line(""" 
        bcftools concat -o {output.vcf} {input.vcf} 
            && 
        bgzip -c {output.vcf} >> {output.vcf_gz} 
        """)


# symlink vcf file to maf directory to run vcf2maf
rule _varscan_symlink_maf:
    input:
        vcf = str(rules._varscan_combine_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _varscan_output_vcf:
    input:
        vcf = str(rules._varscan_combine_vcf.output.vcf_gz)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.{vcf_name}.vcf.gz"
    run:
        op.relative_symlink(input.vcf, output.vcf)


rule _varscan_output_maf:
    input:
        vcf = str(rules._varscan_combine_vcf.output.vcf), # ensure vcf is not deleted before maf is created
        maf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.maf"
    output:
        maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.{vcf_name}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)


def _varscan_get_output(wildcards):
    return expand([
        str(rules._varscan_output_vcf.output.vcf),
        str(rules._varscan_output_maf.output.maf)
        ],
        vcf_name = ["snp", "indel"], **wildcards)


rule _varscan_dispatch:
    input:
        _varscan_get_output
    output:
        dispatched = touch(CFG["dirs"]["outputs"] + "dispatched/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.dispatched")


# Generates the target sentinels for each run, which generate the symlinks
rule _varscan_all:
    input:
        expand(
            [
                str(rules._varscan_dispatch.output.dispatched),
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
