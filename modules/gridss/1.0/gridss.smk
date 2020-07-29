#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["gridss"]`
CFG = op.setup_module(
    name = "gridss",
    version = "1.0",
    subdirectories = ["inputs", "gridss", "viral_annotation", "somatic_filter", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gridss_input_bam,
    _gridss_input_references,
    _gridss_somatic_filter,
    _gridss_output_vcf,
    _gridss_all

wildcard_constraints: 
    genome_build = "|".join(CFG["switches"]["pon"].keys())


##### RULES #####

rule _gridss_input_references: 
    input: 
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_bwa_prefix = reference_files("genomes/{genome_build}/bwa_index/bwa_index/bwa-0.7.17/genome.fa")
    output: 
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome_fa/genome.fa"
    shell: 
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fa}.fai {input.genome_fa}.fai &&
        ln -s {input.genome_bwa_prefix} `dirname {output.genome_fa}`
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gridss_input_bam:
    input:
        sample_bam = CFG["inputs"]["sample_bam"], 
        sample_bai = CFG["inputs"]["sample_bai"] 
    output:
        sample_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam", 
        sample_bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai" 
    run:
        op.relative_symlink(input.sample_bam, output.sample_bam)
        op.relative_symlink(input.sample_bai, output.sample_bai)


# Run GRIDSS paired or unpaired mode
rule _gridss_paired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_id}.bam",
        fasta = rules._gridss_input_references.output.genome_fa, 
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.bed"), 
        repeatmasker = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        vcf = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/gridss.vcf.gz",
        assembly = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/assembly.bam"
    log: CFG["logs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--matched/gridss.log"
    params:
        opts = CFG["options"]["gridss"], 
        steps = "preprocess,assemble,call" 
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["gridss"]
    resources:
        mem_mb = CFG["mem_mb"]["gridss"], 
        bam = 1
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --output {output.vcf}
        --workingdir `dirname {output.vcf}`
        --assembly {output.assembly}
        --blacklist {input.blacklist}
        --repeatmaskerbed {input.repeatmasker}
        --threads {threads}
        --jvmheap {resources.mem_mb}m
        --labels "{wildcards.normal_id},{wildcards.tumour_id}"
        --steps {params.steps}
        {params.opts}
        {input.normal_bam} 
        {input.tumour_bam} 
        2>&1 | tee -a {log}
        """)

rule _gridss_unpaired:
    input:
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        fasta = rules._gridss_input_references.output.genome_fa, 
        blacklist = reference_files("genomes/{genome_build}/encode/encode-blacklist.bed"), 
        repeatmasker = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        vcf = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/gridss.vcf.gz",
        assembly = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/assembly.bam"
    log: CFG["logs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--None--no_normal/gridss.log"
    params:
        opts = CFG["options"]["gridss"], 
        steps = "preprocess,assemble,call"
    conda:
        CFG["conda_envs"]["gridss"]
    threads:
        CFG["threads"]["gridss"]
    resources:
        mem_mb = CFG["mem_mb"]["gridss"], 
        bam = 1
    shell:
        op.as_one_line("""
        gridss
        --reference {input.fasta}
        --output {output.vcf}
        --workingdir `dirname {output.vcf}`
        --assembly {output.assembly}
        --blacklist {input.blacklist}
        --repeatmaskerbed {input.repeatmasker}
        --threads {threads}
        --jvmheap {resources.mem_mb}m
        --labels "{wildcards.tumour_id}"
        --steps {params.steps}
        {params.opts}
        {input.tumour_bam} 
        2>&1 | tee -a {log}
        """)


# Perform viral annotation of the output VCFs
rule _gridss_viral_annotation: 
    input: 
        vcf = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss.vcf.gz", 
        tmpdir = CFG["dirs"]["gridss"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.bam.gridss.working"
    output: 
        vcf = CFG["dirs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation.vcf.gz", 
    log: CFG["logs"]["viral_annotation"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_viral_annotation.log"
    params: 
        viral_ref = CFG["references"]["viral_fa"],
        gridss_jar = "$(readlink -e gridss).jar"
    conda: 
        CFG["conda_envs"]["gridss"]
    threads: 
        CFG["threads"]["viral_annotation"]
    resources: 
        mem_mb = CFG["mem_mb"]["viral_annotation"]
    shell: 
        op.as_one_line("""
        java -Xmx{resources.mem_mb}m 
                -cp {params.gridss_jar} gridss.AnnotateInsertedSequence 
				REFERENCE_SEQUENCE={params.viral_ref} 
				INPUT={input.vcf} 
				OUTPUT={output.vcf} 
				WORKER_THREADS={threads} 
                2>&1 | tee -a {log}
        """)
    
rule _gridss_somatic_filter: 
    input: 
        vcf = rules._gridss_viral_annotation.output.vcf
    output: 
        somatic_vcf = CFG["dirs"]["somatic_filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic.vcf.gz", 
        full_vcf = CFG["dirs"]["somatic_filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_full.vcf.gz"
    log: log = CFG["logs"]["somatic_filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/gridss_somatic_filter.log"
    params:
        pondir = op.switch_on_wildcard("genome_build", CFG["switches"]["pon"]),
        bsgenome = lambda w: {
            "grch37": "hg19", 
            "hs37d5": "hg19", 
            "hg38": "hg38"}[w.genome_build], 
        scriptdir = "$(dirname $(readlink -e $(which gridss)))"
    shell: 
        op.as_one_line(""" 
        gridss_somatic_filter 
        -p {params.pondir} 
        -i {input.vcf}
        -o {output.somatic_vcf}
        -f {output.full_vcf}
        -r BSgenome.Hsapiens.UCSC.{params.bsgenome}
        -s {params.scriptdir}
        --gc 
        2>&1 | tee -a {log}
        """)
    

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gridss_output_vcf:
    input:
        somatic = str(rules._gridss_somatic_filter.output.somatic_vcf), 
        full = str(rules._gridss_somatic_filter.output.full_vcf)
    output:
        somatic = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic.vcf.gz", 
        full = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.gridss_somatic_full.vcf.gz"
    run:
        op.relative_symlink(input.somatic, output.somatic)
        op.relative_symlink(input.full, output.full)


# Generates the target sentinels for each run, which generate the symlinks
rule _gridss_all:
    input:
        expand(
            [
                str(rules._gridss_output_vcf.output.somatic),
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
