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
    subdirectories = ["inputs", "varscan", "maf", "outputs"]
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
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



##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _varscan_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bam + ".bai", output.bam + ".bai")

# Pulls in list of chromosomes for the genome builds
checkpoint _varscan_input_chroms:
    input:
        txt = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        txt = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.txt, output.txt)

#generate mpileups for tumour and normal bams separately. 
#If we parallelize this by chromosome we will need 2 * 2 threads per chromosome but this should be a lot more efficient
rule _varscan_bam2mpu:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        mpu = pipe(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{sample_id}.{chrom}.mpileup")
    log:
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{sample_id}.bam2mpu.{chrom}.stderr.log"
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
        tumour_mpu = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}.{chrom}.mpileup",
        normal_mpu = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{normal_id}.{chrom}.mpileup"
    output:
        snp = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.snp.vcf"),
        indel = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.indel.vcf")
    wildcard_constraints:
        pair_status = "matched|unmatched"
    log:
        stdout = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stdout.log",
        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.varscan_somatic.stderr.log"
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
#added because varscan seemed to be exiting without error  

rule _varscan_reheader_vcf:
    input:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.{vcf_name}.vcf",
        #snp = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.snp.vcf",
        #indel = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.indel.vcf",
        header = CFG["vcf_header"]  #need to make this work for different genome builds
    output:
        vcf = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.{vcf_name}.vcf.gz"),
        #snp = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.snp.vcf.gz"),
        #indel = temp(CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{chrom}.indel.vcf.gz")
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools reheader -h {input.header} {input.vcf} | bcftools view -l 1 -o {output.vcf}
        """)

#bcftools reheader -h {input.header} {input.snp} | bcftools view -l 1 -o {output.snp} && bcftools reheader -h {input.header} {input.indel} | bcftools view -l 1 -o {output.indel} 

#rule _varscan_unpaired:
#    input:
#        tumour_mpu = CFG["dirs"]["mpileup"] + "{seq_type}--{genome_build}/{tumour_id}.mpileup",
#    output:
#        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf"
#    wildcard_constraints:
#        pair_status = "no_normal"
#    log:
#        stderr = CFG["logs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/varscan_{vcf_name}.stderr.log"
#    params:
#        opts = op.switch_on_wildcard("seq_type", CFG["options"]["unpaired"]),
#        cns = op.switch_on_wildcard("seq_type", {"cns": "--variants", "indel": "", "snp": ""})
#    conda:
#        CFG["conda_envs"]["varscan"]
#    threads:
#        CFG["threads"]["unpaired"]
#    resources:
#        mem_mb = CFG["mem_mb"]["unpaired"]
#    shell:
#        op.as_one_line("""
#        varscan mpileup2{wildcards.vcf_name} 
#        {input.tumour_mpu} 
#        {params.opts}
#        {params.cns}
#        > {output.vcf} 2> {log.stderr}
#        """)


def _varscan_request_chrom_vcf(wildcards):
    CFG = config["lcr-modules"]["varscan"]
    with open(checkpoints._varscan_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    print(mains_chroms)
    #the line below was causing and error and I don't know why
    vcf_files = expand(rules._varscan_reheader_vcf.output.vcf,
        chrom = mains_chroms, **wildcards
    )
    print(vcf_files)
    return vcf_files

"""
def _varscan_request_chrom_vcf_indel(wildcards):
    CFG = config["lcr-modules"]["varscan"]
    with open(checkpoints._varscan_input_chroms.get(**wildcards).output.txt) as f:
        mains_chroms = f.read().rstrip("\n").split("\n")
    print(mains_chroms)
    #the line below was causing and error and I don't know why
    #vcf_files = expand(
    #    CFG["dirs"]["varscan"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{chrom}.{vcf_name}.vcf", 
    #    chrom=mains_chroms, vcf_name=['snp','indel']
    #)
    vcf_files = expand("results/varscan-1.0/01-varscan/{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{normal_id}}--{{pair_status}}/{chrom}.indel.vcf.gz",
     chrom=mains_chroms)
    print(vcf_files)
    return vcf_files
"""

# I generalize this rule for any vcf_type
rule _varscan_combine_vcf:
    input:
        vcf = _varscan_request_chrom_vcf
    output:
        vcf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pass.somatic.{vcf_name}.vcf.gz",
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        " bcftools concat -o {output.vcf} {input.vcf}"

#currently disabled
rule _varscan_symlink_maf:
    input:
        vcf = rules._varscan_combine_vcf.output.vcf
        #CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz"
    output:
        vcf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.vcf.gz"
    run:
        op.relative_symlink(input.vcf, output.vcf)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _varscan_output_vcf:
    input:
        vcf = rules._varscan_combine_vcf.output.vcf
        #snp_vcf = rules._varscan_combine_vcf.output.snp,
        #indel_vcf = rules._varscan_combine_vcf.output.indel
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.{vcf_name}.vcf.gz"
        #indel_vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.indel.vcf.gz",
        #snp_vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.snp.vcf.gz"
    run:
        op.relative_symlink(input.vcf, output.vcf)

#currently not working. Add to inputs at end?
rule _varscan_output_maf:
    input:
        maf = CFG["dirs"]["maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{vcf_name}.maf"
        #snp_maf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pass.somatic.snp.maf",
        #indel_maf = CFG["dirs"]["varscan"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/pass.somatic.indel.maf"
    output:
        maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}-pass.somatic.{vcf_name}.maf"
        #snp_maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/snp/{tumour_id}--{normal_id}--{pair_status}.pass.somatic.snp.maf",
        #indel_maf = CFG["dirs"]["outputs"] + "maf/{seq_type}--{genome_build}/snp/{tumour_id}--{normal_id}--{pair_status}.pass.somatic.indel.maf"
    run:
        op.relative_symlink(input.maf, output.maf)


def _varscan_get_output(wildcards):
    if wildcards.pair_status == "no_normal":
        #return expand([rules._varscan_combine_vcf.output.vcf, rules._varscan_output_maf.output.maf], vcf_name = ["indel", "snp", "cns"], **wildcards)
        raise ValueError("wildcards.pair_status = non_normal is currently unsupported, please run in matched or unmatched mode")
        #return expand([rules._varscan_combine_vcf.output.snp, rules._varscan_combine_vcf.output.indel], **wildcards)
    else: 
        return expand([
            rules._varscan_output_vcf.output.vcf,rules._varscan_output_maf.output.maf
            ], vcf_name = ["snp", "indel"], **wildcards)
        #return expand([rules._varscan_combine_vcf.output.{vcf_name}, rules._varscan_output_maf.output.maf], vcf_name = ["indel", "snp"], **wildcards)

#I'm not sure how to change this to not need to expand the chromosome names. Is this to deal with checkpoints? 
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
                rules._varscan_dispatch.output.dispatched,
                #rules._varscan_output.output.snp_vcf,
                #rules._varscan_output_maf.output.snp_maf,
                #rules._varscan_output_maf.output.indel_maf
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            #chrom = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
