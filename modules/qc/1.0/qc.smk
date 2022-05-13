#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostia Dreval
# Module Author:    Kostia Dreval
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["qc"]`
CFG = op.setup_module(
    name = "qc",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "samtools", "gatk", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _qc_input_bam,
    _qc_step_2,
    _qc_output_tsv,
    _qc_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
# TODO: If applicable, add an input rule for each input file used by the module
# TODO: If applicable, create second symlink to .crai file in the input function, to accomplish cram support
rule _qc_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)

# symlink the reference files to ensure all index/dictionaries are available for GATK tools
rule _qc_input_references:
    input:
        genome_fa = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        genome_fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        genome_dict = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        genome_fa = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa",
        genome_fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai",
        genome_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.dict"
    shell:
        op.as_one_line("""
        ln -s {input.genome_fa} {output.genome_fa} &&
        ln -s {input.genome_fai} {output.genome_fai} &&
        ln -s {input.genome_dict} {output.genome_dict}
        """)


# Collect samtools stats
rule _qc_samtools_stat:
    input:
        bam = str(rules._qc_input_bam.output.bam)
    output:
        samtools_stat = CFG["dirs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.stat"
    log:
        stderr = CFG["logs"]["samtools"] + "{seq_type}--{genome_build}/{sample_id}.run_samtools.stderr.log"
    params:
        opts = CFG["options"]["samtools_stat"]
    conda:
        CFG["conda_envs"]["samtools"]
    threads:
        CFG["threads"]["samtools_stat"]
    resources:
        **CFG["resources"]["samtools_stat"]
    shell:
        op.as_one_line("""
        samtools stat
        {params.opts}
        --threads {threads}
        {input.bam}
        >> {output.flagstat}
        2>> {log.stderr}
        """)


rule _qc_gatk_basequality:
    input:
        bam = str(rules._qc_input_bam.output.bam)
    output:
        gatk_basequal = CFG["dirs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.txt",
        gatk_basequal_chart = CFG["dirs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.pdf"
    log:
        stdout = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stdout.log",
        stderr = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stderr.log"
    params:
        opts = CFG["options"]["QualityScoreDistribution"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    threads:
        CFG["threads"]["QualityScoreDistribution"]
    resources:
        **CFG["resources"]["QualityScoreDistribution"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        gatk QualityScoreDistribution --spark-master local[{threads}]
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        {params.opts}
        -I {input.bam}
        -O {output.gatk_basequal}
        -CHART {output.gatk_basequal_chart}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)

rule _qc_gatk_wgs:
    input:
        bam = str(rules._qc_input_bam.output.bam),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        gatk_wgs = CFG["dirs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectWgsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stderr.log"
    wildcard_constraints:
        seq_type = "genome"
    params:
        opts = CFG["options"]["CollectWgsMetrics"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    threads:
        CFG["threads"]["CollectMetrics"]
    resources:
        **CFG["resources"]["CollectMetrics"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        gatk CollectWgsMetrics --spark-master local[{threads}]
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        {params.opts}
        -I {input.bam}
        -O {output.gatk_basequal}
        -R {input.fasta}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


rule _qc_gatk_wes:
    input:
        bam = str(rules._qc_input_bam.output.bam),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        gatk_wes = CFG["dirs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectHsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stderr.log"
    wildcard_constraints:
        seq_type = "capture"
    params:
        opts = CFG["options"]["CollectHsMetrics"],
        jvmheap = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatkR"]
    threads:
        CFG["threads"]["CollectMetrics"]
    resources:
        **CFG["resources"]["CollectMetrics"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        gatk CollectHsMetrics --spark-master local[{threads}]
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        {params.opts}
        -I {input.bam}
        -O {output.gatk_basequal}
        -R {input.fasta}
        -BI
        -TI
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)



# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
rule _qc_step_2:
    input:
        tsv = str(rules._qc_step_1.output.tsv)
    output:
        tsv = CFG["dirs"]["qc"] + "{seq_type}--{genome_build}/{sample_id}/output.filt.tsv"
    log:
        stderr = CFG["logs"]["qc"] + "{seq_type}--{genome_build}/{sample_id}/step_2.stderr.log"
    params:
        opts = CFG["options"]["step_2"]
    shell:
        "grep {params.opts} {input.tsv} > {output.tsv} 2> {log.stderr}"


# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _qc_output_tsv:
    input:
        tsv = str(rules._qc_step_2.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{seq_type}--{genome_build}/{sample_id}.output.filt.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _qc_all:
    input:
        expand(
            [
                str(rules._qc_output_tsv.output.tsv),
                # TODO: If applicable, add other output rules here
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
