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
    subdirectories = ["inputs", "samtools", "gatk", "aggregated_metrics", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _qc_input_bam,
    _qc_input_references,
    _qc_output_tsv,
    _qc_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _qc_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bam + ".bai", output.bai)
        op.relative_symlink(input.bam + ".bai", output.crai)

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
        samtools stats
        {params.opts}
        --threads {threads}
        {input.bam}
        >> {output.samtools_stat}
        2>> {log.stderr}
        """)

# Collecting GATK base quality metrics
rule _qc_gatk_basequality:
    input:
        bam = str(rules._qc_input_bam.output.bam),
        fasta = str(rules._qc_input_references.output.genome_fa)
    output:
        gatk_basequal = CFG["dirs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.txt",
        gatk_basequal_chart = CFG["dirs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.{genome_build}.QualityScoreDistribution.pdf"
    log:
        stdout = CFG["logs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stdout.log",
        stderr = CFG["logs"]["gatk"] + "QualityScoreDistribution/{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.stderr.log"
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
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        QualityScoreDistribution
        {params.opts}
        -I {input.bam}
        -O {output.gatk_basequal}
        -CHART {output.gatk_basequal_chart}
        -R {input.fasta}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


# Collect GATK QC metrics (separately for genomes/exomes)
rule _qc_gatk_wgs:
    input:
        bam = str(rules._qc_input_bam.output.bam),
        fasta = str(rules._qc_input_references.output.genome_fa),
        samtools_stats = str(rules._qc_samtools_stat.output.samtools_stat)
    output:
        gatk_wgs = CFG["dirs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectWgsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectWgsMetrics.stderr.log"
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
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        CollectWgsMetrics
        {params.opts}
        -I {input.bam}
        -O {output.gatk_wgs}
        -R {input.fasta}
        --READ_LENGTH $(grep ^SN {input.samtools_stats} | cut -f 2- | grep "average length:" | cut -f 2)
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


def _qc_get_baits(wildcards):
    CFG = config["lcr-modules"]["qc"]
    # simplify config options by re-using the capture space bertween e.g grch37 and hs37d5
    if "38" in str(wildcards.genome_build):
        this_genome_build = "hg38"
    else:
        this_genome_build = "grch37"

    if "baits_regions" in  CFG:
        # If user specified path to bed file in the sample table, use it from dictionary in config
        if wildcards.baits_regions in list(CFG["baits_regions"][this_genome_build].keys()):
            these_baits = str(wildcards.baits_regions)
        else:
            print(
                f"WARNING: the baits regions were specified in the sample table, but were not found in config "
                f"for {wildcards.genome_build} ({this_genome_build}) and {wildcards.seq_type} combination. Using the default space ..."
            )
            these_baits = "_default"
    else:
        # use default
        these_baits = "_default"

    return(str(CFG["baits_regions"][this_genome_build][these_baits]))


# Subset contigs in bed to those present in the reference and sort the output
# This is needed to pass GATK validation
rule _qc_sort_baits:
    input:
        baits = _qc_get_baits,
        fai = str(rules._qc_input_references.output.genome_fai)
    output:
        baits = CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.bed",
        inermediate_baits = temp(CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.INTEMEDIATE.bed")
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        if [ -e {input.baits} ]; then
            cat {input.baits} > {output.inermediate_baits};
        else
            curl -L {input.baits} > {output.inermediate_baits};
        fi
            &&
        cut -f 1-2  {input.fai} |
        perl -ane 'print "$F[0]\\t0\\t$F[1]\\n"' |
        bedtools intersect -wa -a {output.inermediate_baits} -b stdin |
        bedtools sort |
        bedtools merge >
        {output.baits}
        """)

# Create interval list
rule _qc_baits_to_interval_list:
    input:
        baits = str(rules._qc_sort_baits.output.baits),
        genome_dict = str(rules._qc_input_references.output.genome_dict)
    output:
        interval_list = CFG["dirs"]["inputs"] + "references/{genome_build}/{baits_regions}.interval_list"
    conda:
        CFG["conda_envs"]["gatkR"]
    shell:
        op.as_one_line("""
        gatk BedToIntervalList
        --INPUT {input.baits}
        -SD {input.genome_dict}
        -O {output.interval_list}
        """)

# get the proper interval list corresponding to each sample, if provided
def _qc_get_intervals(wildcards):
    CFG = config["lcr-modules"]["qc"]
    samples = CFG["samples"]
    this_sample = samples.loc[(samples['sample_id'] == wildcards.sample_id) &
            (samples['genome_build'] == wildcards.genome_build) &
            (samples['seq_type'] == wildcards.seq_type)]

    # ensure the sample is uniquely mapped
    if len(this_sample) != 1:
        raise AssertionError("Found %s matches when examining the sample table for \'%s\' \'%s\' \'%s\'" % (len(sample), sample_id, genome_build, seq_type))

    # simplify config options by re-using the capture space bertween e.g grch37 and hs37d5
    if "38" in str(wildcards.genome_build):
        this_genome_build = "hg38"
    else:
        this_genome_build = "grch37"

    if "baits_regions" in  this_sample:
        # If user specified path to bed file in the sample table, use it from dictionary in config
        if str(this_sample.iloc[0]['baits_regions']) in list(CFG["baits_regions"][this_genome_build].keys()):
            these_baits = str(this_sample.iloc[0]['baits_regions'])
        else:
            print(
                f"WARNING: the baits regions were specified in the sample table, but were not found in config "
                f"for {wildcards.genome_build} ({this_genome_build}) and {wildcards.seq_type} combination. Using the default space ..."
            )
            these_baits = "_default"
    else:
        # use default
        these_baits = "_default"
    return(expand(str(rules._qc_baits_to_interval_list.output.interval_list), baits_regions=these_baits,allow_missing=True))

# Collect metrics on WES samples
rule _qc_gatk_wes:
    input:
        bam = str(rules._qc_input_bam.output.bam),
        fasta = str(rules._qc_input_references.output.genome_fa),
        intervals = _qc_get_intervals
    output:
        gatk_wes = CFG["dirs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.{genome_build}.CollectHsMetrics.txt"
    log:
        stdout = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stdout.log",
        stderr = CFG["logs"]["gatk"] + "CollectMetrics/{seq_type}--{genome_build}/{sample_id}.CollectHsMetrics.stderr.log"
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
        gatk
        --java-options "-Xmx{params.jvmheap}m -XX:ConcGCThreads=1"
        CollectHsMetrics
        {params.opts}
        -I {input.bam}
        -O {output.gatk_wes}
        -R {input.fasta}
        -BI {input.intervals}
        -TI {input.intervals}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.sample_id} on $(hostname) at $(date)" >> {log.stdout};
        """)


def _qc_get_stats(wildcards):
    if wildcards.seq_type in ["capture"]:
        output = str(rules._qc_gatk_wes.output.gatk_wes)
    else:
        output = str(rules._qc_gatk_wgs.output.gatk_wgs)
    return output


# Collect required metrix into a tidy table
rule _qc_collect_metrics:
    input:
        stat = str(rules._qc_samtools_stat.output.samtools_stat),
        base_qual = str(rules._qc_gatk_basequality.output.gatk_basequal),
        metrics = _qc_get_stats
    output:
        stat = CFG["dirs"]["aggregated_metrics"] + "{seq_type}--{genome_build}/{sample_id}.metrix.tsv"
    conda:
        CFG["conda_envs"]["gatkR"]
    threads:
        CFG["threads"]["collect"]
    resources:
        **CFG["resources"]["collect"]
    script:
        "src/R/aggregate_metrics.R"


# Combine all individual metrics into one single file
rule _qc_merge_metrics:
    input:
        expand(
            [
                str(rules._qc_collect_metrics.output.stat)
            ],
            zip,
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"],
            allow_missing = True)
    output:
        CFG["dirs"]["aggregated_metrics"] + "{seq_type}/master_qc_metrics_output.tsv"
    shell:
        op.as_one_line("""
        head -n1 {input[0]} > {output};
        tail -n+2 -q {input} >> {output};
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _qc_symlink_output:
    input:
        str(rules._qc_merge_metrics.output)
    output:
        CFG["dirs"]["outputs"] + "{seq_type}.qc_metrics.tsv",
    run:
        op.relative_symlink(input, output, in_module= True)

rule _qc_output_tsv:
    input:
        stat = str(rules._qc_samtools_stat.output.samtools_stat),
        base_qual = str(rules._qc_gatk_basequality.output.gatk_basequal),
        metrics = _qc_get_stats
    output:
        stat = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.stat.tsv",
        base_qual = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.QualityScoreDistribution.tsv",
        metrics = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{sample_id}.CollectMetrics.tsv",
    run:
        op.relative_symlink(input.stat, output.stat, in_module= True)
        op.relative_symlink(input.base_qual, output.base_qual, in_module= True)
        op.relative_symlink(input.metrics, output.metrics, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _qc_all:
    input:
        expand(
            [
                str(rules._qc_output_tsv.output.stat),
                str(rules._qc_output_tsv.output.base_qual),
                str(rules._qc_output_tsv.output.metrics)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["samples"]["seq_type"],
            genome_build=CFG["samples"]["genome_build"],
            sample_id=CFG["samples"]["sample_id"]),
        expand(
            [
                str(rules._qc_symlink_output.output)
            ],
            seq_type=list(CFG["samples"]["seq_type"].unique()))


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
