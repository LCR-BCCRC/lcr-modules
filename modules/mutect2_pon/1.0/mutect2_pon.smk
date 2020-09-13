#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Prasath Pararajalingam
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mutect2_pon"]`
CFG = op.setup_module(
    name = "mutect2_pon",
    version = "1.0",
    subdirectories = ["inputs", "mutect2", "genomicsdb", "panel_of_normals", "outputs"]
)

# Create a collapsed sample table to link sample_id to genome_build and seq_type
SAMPLES_COLLAPSED = CFG["samples"].groupby(["genome_build", "seq_type"])["sample_id"].apply(list)
SAMPLES_COLLAPSED = SAMPLES_COLLAPSED.reset_index()

# Define rules to be run locally when using a compute cluster
localrules:
    _mutect2_pon_input_bam,
    _mutect2_pon_input_chrs,
    _mutect2_pon_get_sm,
    _mutect2_pon_output_vcf,
    _mutect2_pon_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _mutect2_pon_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)


# Symlink chromosomes used for parallelization
checkpoint _mutect2_pon_input_chrs:
    input:
        chrs = reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.relative_symlink(input.chrs, output.chrs)


# Retrieves from SM tag from BAM and writes to file
rule _mutect2_pon_get_sm:
    input:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam"
    output:
        sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}_sm.txt",
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}_mutect2_pon_get_sm.stderr.log"
    conda:
        CFG["conda_envs"]["samtools"]
    shell:
        "samtools view -H {input.bam} | grep '^@RG' | "
        r"sed 's/.*SM:\([^\t]*\).*/\1/g'"" | uniq > {output.sm} 2> {log.stderr}"

# Launces Mutect2 in no normal mode
checkpoint _mutect2_pon_run:
    input:
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        normal_sm = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}_sm.txt"
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/chromosomes/{chrom}.output.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/chromosomes/{chrom}.output.vcf.gz.tbi"),
        stat = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/chromosomes/{chrom}.output.vcf.gz.stats")
    log:
        stdout = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.mutect2_run.stdout.log",
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/{chrom}.mutect2_run.stderr.log"
    resources:
        **CFG["resources"]["mutect2_run"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2_run"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_run"]
    wildcard_constraints: 
        pair_status = "no_normal"
    shell:
        op.as_one_line("""
        gatk Mutect2 --java-options "-Xmx{params.mem_mb}m" {params.opts} -I {input.normal_bam} -R {input.fasta} 
        -O {output.vcf} --germline-resource {input.gnomad} -L {wildcards.chrom}
        > {log.stdout} 2> {log.stderr}
        """)


def _mutect2_pon_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["mutect2_pon"]
    chrs = checkpoints._mutect2_pon_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    vcfs = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}/chromosomes/{chrom}.output.vcf.gz",
        chrom = chrs
    )
    return(vcfs)


def _mutect2_pon_get_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["mutect2_pon"]
    chrs = checkpoints._mutect2_pon_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    tbis = expand(
        CFG["dirs"]["mutect2"] + "{{seq_type}}--{{genome_build}}/{{sample_id}}/chromosomes/{chrom}.output.vcf.gz.tbi",
        chrom = chrs
    )
    return(tbis)


# Merge chromosome mutect2 VCFs from the same sample
rule _mutect2_pon_merge_vcfs:
    input:
        vcf = _mutect2_pon_get_chr_vcfs,
        tbi = _mutect2_pon_get_chr_tbis
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/mutect2_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["mutect2_merge_vcfs"]
    resources:
        **CFG["resources"]["mutect2_merge_vcfs"]
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -O z -o {output.vcf} 2> {log.stderr} 
            &&
        bcftools index -t --threads {threads} {output.vcf} 2>> {log.stderr}
        """)

def _mutect2_pon_gather_sms(wildcards):
    CFG = config["lcr-modules"]["mutect2_pon"]
    sample_ids = CFG["samples"][(CFG["samples"].genome_build == wildcards.genome_build) & (CFG["samples"].seq_type == wildcards.seq_type)]["sample_id"]
    sms = expand(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}_sm.txt", 
                    sample_id = sample_ids, 
                    seq_type = wildcards.seq_type, 
                    genome_build = wildcards.genome_build
    )
    return sms 
    
def _mutect2_pon_gather_vcfs(wildcards): 
    CFG = config["lcr-modules"]["mutect2_pon"]
    sample_ids = CFG["samples"][(CFG["samples"].genome_build == wildcards.genome_build) & (CFG["samples"].seq_type == wildcards.seq_type)]["sample_id"]
    vcfs = expand(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz", 
                    sample_id = sample_ids, 
                    seq_type = wildcards.seq_type, 
                    genome_build = wildcards.genome_build
    )
    return vcfs

def _mutect2_pon_gather_tbis(wildcards): 
    CFG = config["lcr-modules"]["mutect2_pon"]
    sample_ids = CFG["samples"][(CFG["samples"].genome_build == wildcards.genome_build) & (CFG["samples"].seq_type == wildcards.seq_type)]["sample_id"]
    tbis = expand(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}/output.vcf.gz.tbi", 
                    sample_id = sample_ids, 
                    seq_type = wildcards.seq_type, 
                    genome_build = wildcards.genome_build
    )
    return tbis
    

# Generate a table of normals with their SM tag
rule _mutect2_pon_vcf_table: 
    input: 
        vcfs = _mutect2_pon_gather_vcfs, 
        tbis = _mutect2_pon_gather_tbis,
        sms = _mutect2_pon_gather_sms   
    output: 
        table = CFG["dirs"]["genomicsdb"] + "{seq_type}--{genome_build}/normal_vcfs.tsv"
    conda: 
        CFG["conda_envs"]["coreutils"]
    shell: 
        op.as_one_line("""
        pr -mts <(cat {input.sms}) <(ls {input.vcfs}) > {output.table}
        """)


# Build genomicsDB from all normal VCFs
rule _mutect2_pon_genomicsdb: 
    input:
        vcfs = _mutect2_pon_gather_vcfs, 
        tbis = _mutect2_pon_gather_tbis,
        table = str(rules._mutect2_pon_vcf_table.output.table), 
        chrs = str(rules._mutect2_pon_input_chrs.output.chrs)
    output: 
        intlist = CFG["dirs"]["genomicsdb"] + "{seq_type}--{genome_build}/interval_list.list",
        genomicsdb = CFG["dirs"]["genomicsdb"] + "{seq_type}--{genome_build}/pon_db"
    log: 
        CFG["logs"]["genomicsdb"] + "{seq_type}--{genome_build}/genomicsdb.log" 
    resources: 
        **CFG["resources"]["genomicsdb"]
    params: 
        mem_mb = lambda wildcard, resources: int(resources.mem_mb * 0.8)
    conda: 
        CFG["conda_envs"]["gatk"]
    threads: 
        CFG["threads"]["genomicsdb"]
    shell: 
        op.as_one_line("""
        cp -f {input.chrs} {output.intlist} &&
        gatk --java-options "-Xmx{params.mem_mb}m -Xms{params.mem_mb}m" GenomicsDBImport 
            --genomicsdb-workspace-path $PWD/{output.genomicsdb} 
            --sample-name-map {input.table} 
            --reader-threads {threads}
            -L {output.intlist}
        """)

# Create panel of normals vcf from genomicsDB
rule _mutect2_pon_create_pon: 
    input: 
        genomicsdb = rules._mutect2_pon_genomicsdb.output.genomicsdb, 
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output: 
        vcf = CFG["dirs"]["panel_of_normals"] + "{seq_type}--{genome_build}/pon.vcf.gz", 
        tbi = CFG["dirs"]["panel_of_normals"] + "{seq_type}--{genome_build}/pon.vcf.gz.tbi"
    resources: 
        **CFG["resources"]["create_pon"]
    params: 
        mem_mb = lambda wildcard, resources: int(resources.mem_mb * 0.8)
    conda: 
        CFG["conda_envs"]["gatk"]
    threads: 
        CFG["threads"]["create_pon"]
    shell: 
        op.as_one_line("""
            gatk CreateSomaticPanelOfNormals --java-options "-Xmx{params.mem_mb}m" 
                -R {input.fasta} 
                --germline-resource {input.gnomad} 
                -V gendb:/$PWD/{input.genomicsdb} 
                -O {output.vcf}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _mutect2_pon_output_vcf:
    input:
        vcf = str(rules._mutect2_pon_create_pon.output.vcf),
        tbi = str(rules._mutect2_pon_create_pon.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}.pon.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}.pon.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf)
        op.relative_symlink(input.tbi, output.tbi)


# Generates the target sentinels for each run, which generate the symlinks
rule _mutect2_pon_all:
    input:
        expand(
            [
                str(rules._mutect2_pon_output_vcf.output.vcf),
                str(rules._mutect2_pon_output_vcf.output.tbi)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=SAMPLES_COLLAPSED["seq_type"],
            genome_build=SAMPLES_COLLAPSED["genome_build"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
