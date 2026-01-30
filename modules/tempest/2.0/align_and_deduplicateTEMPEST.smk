#!/usr/bin/env snakemake

"""
TEMPEST cfDNA UMI Duplex Consensus Pipeline v2

Key changes from v1:
  - FastqToBam replaces fastp + AnnotateBamWithUmis (eliminates memory bottleneck)
  - ZipperBams replaces Picard MergeBamAlignment
  - Combined align → sort → group → consensus → filter in single step
  - Tag removal via samtools view -x (replaces sanitize_bam.py)
  - Streaming throughout with --async-io and --compression optimization
  - Dynamic memory allocation with retry scaling

Pipeline structure:
  Step 1: FASTQ → Unmapped BAM (with UMIs)
  Step 2: Unmapped BAM → Consensus Unmapped BAM (align, sort)
  Step 3: Group by UMI
  Step 4: Call and Filter Consensus
  Step 3: Consensus Unmapped BAM → Final BAM (realign, remove tags, sort, index)
"""

import os
import sys
import pandas as pd

sys.path.append(os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "modules/tempest/2.0/"))

BAM_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "tempest/2.0/", "bam_pipeline")
UTILSDIR = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "modules/tempest/2.0/utils")
UMI_SAMPLESHEET = config["lcr-modules"]["_shared"]["samples"].copy()


# =============================================================================
# Helper Functions
# =============================================================================

def get_capture_regions(wildcards):
    """Get the capture regions file path for a sample from the sample sheet and config."""
    capture_space = UMI_SAMPLESHEET.loc[UMI_SAMPLESHEET["sample_id"] == wildcards.sample, "capture_space"]
    if capture_space.empty:
        raise ValueError(f"Capture space not found for {wildcards.sample}")
    capture_space_value = capture_space.values[0]
    if capture_space_value not in config["lcr-modules"]["_shared"]["captureregionsil"]:
        raise ValueError(f"Capture space '{capture_space_value}' not found in config")
    return config["lcr-modules"]["_shared"]["captureregionsil"][capture_space_value]


# =============================================================================
# Dynamic Memory Functions
# =============================================================================

def fastq_to_bam_mem(wildcards, input, attempt):
    """
    Memory for FastqToBam step.
    Streams through FASTQs, so memory needs are minimal and stable.
    """
    if attempt == 1:
        return 4000
    elif attempt == 2:
        return 6000
    elif attempt == 3:
        return 10000
    else:
        return 16000

def bwa_mem(wildcards, input, attempt):
    """
    Memory for align + ZipperBams step.
    Scales with input unmapped BAM size.
    """
    if attempt == 1:
        return 15000
    elif attempt == 2:
        return max(input.size_mb, 20000)
    elif attempt == 3:
        return min(input.size_mb * 1.5, 50000)
    else:
        return min(input.size_mb * 2, 80000)

def fgbio_mem(wildcards, input, attempt):
    if attempt == 1:
        return 10000
    elif attempt == 2:
        return max(input.size_mb, 15000)
    elif attempt == 3:
        return min(input.size_mb + 10000, 50000)
    else:
        return min(input.size_mb + 15000, 80000)


def realign_finalize_mem(wildcards, input, attempt):
    """
    Memory for realign → ZipperBams → sort → final BAM step.
    Input is consensus unmapped BAM (typically 10-30% of raw data after collapsing).
    """
    input_size_gb = input.size_mb / 1000
    
    if attempt == 1:
        base_mb = 10000
        per_gb = 150
        cap = 48000
    elif attempt == 2:
        base_mb = 16000
        per_gb = 300
        cap = 64000
    elif attempt == 3:
        base_mb = 24000
        per_gb = 500
        cap = 96000
    else:
        base_mb = 36000
        per_gb = 700
        cap = 128000
    
    calculated = base_mb + (input_size_gb * per_gb)
    return int(min(calculated, cap))


# =============================================================================
# STEP 1: FASTQ -> Unmapped BAM with UMI tags
# =============================================================================

rule fgbio_fastq_to_bam:
    """
    Convert FASTQs to unmapped BAM, extracting UMIs into RX tag.
    
    Replaces: fastp trimming + later AnnotateBamWithUmis
    Memory:   ~4GB (streaming, vs 35-100GB+ for AnnotateBamWithUmis)
    
    Read structures must be specified for both R1 and R2.
    Example: "8M+T 8M+T" = 8bp UMI at start of both reads
    """
    input:
        r1 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r1_path"],
        r2 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r2_path"]
    output:
        bam = temp(os.path.join(BAM_OUTDIR, "01-unmapped", "{sample}.unmapped.bam"))
    params:
        read_structures = config["lcr-modules"]["cfDNA_umi_workflow"]["read_structure"]
    threads: 2
    resources:
        mem_mb = fastq_to_bam_mem
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.fastq_to_bam.log")
    shell:
        """
    fgbio -Xmx{resources.mem_mb}m --compression 1 --async-io FastqToBam \
        --input {input.r1} {input.r2} \
        --read-structures {params.read_structures} \
        --sample {wildcards.sample} \
        --library {wildcards.sample} \
        --output {output.bam} &> {log}
        """

# =============================================================================
# STEP 2: Unmapped BAM -> Consensus Unmapped BAM
#         (align, sort, group, call consensus, filter)
# =============================================================================

rule align_and_zipper:
    """
    Align reads and merge with unmapped BAM tags using ZipperBams.
    """
    input:
        unmapped_bam = rules.fgbio_fastq_to_bam.output.bam,
    output:
        bam = temp(os.path.join(BAM_OUTDIR, "02-aligned", "{sample}.aligned.bam"))
    params:
        bwa_ref = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"],
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["bwa_threads"]
    resources:
        mem_mb = bwa_mem
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.align_and_zipper.log")
    shell:
        """
        samtools fastq {input.unmapped_bam} 2>> {log} \
        | bwa mem -t {threads} -p -K 150000000 -Y {params.bwa_ref} - 2>> {log} \
        | fgbio -Xmx4g --compression 0 --async-io ZipperBams \
            --unmapped {input.unmapped_bam} \
            --ref {params.ref_w_dict} 2>> {log} \
        > {output.bam}
        """

rule group_reads_by_umi:
    """
    Group reads by position + UMI into families.
    """
    input:
        bam = rules.align_and_zipper.output.bam
    output:
        bam = temp(os.path.join(BAM_OUTDIR, "03-grouped", "{sample}.grouped.bam")),
        histogram = os.path.join(BAM_OUTDIR, "03-grouped", "{sample}.family_sizes.txt")
    params:
        grouping_strategy = config["lcr-modules"]["cfDNA_umi_workflow"].get("grouping_strategy", "paired"),
        edits = config["lcr-modules"]["cfDNA_umi_workflow"]["umiedits"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    threads: 2
    resources:
        mem_mb = fgbio_mem
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.group_reads_by_umi.log")
    shell:
        """
          fgbio -Xmx8g --compression 0 --async-io GroupReadsByUmi \
            --input {input.bam} \
            --output {output.bam} \
            --strategy {params.grouping_strategy} \
            --edits {params.edits} \
            --family-size-histogram {output.histogram} 2>> {log}
        """

rule call_and_filter_consensus:
    """
    Generate consensus from families and filter low-quality bases.
    """
    input:
        bam = rules.group_reads_by_umi.output.bam
    output:
        bam = temp(os.path.join(BAM_OUTDIR, "04-consensus", "{sample}.consensus.unmapped.bam"))
    params:
        min_reads = config["lcr-modules"]["cfDNA_umi_workflow"].get("min_reads", "0"),
        min_input_base_quality = config["lcr-modules"]["cfDNA_umi_workflow"].get("min_input_base_quality", 20),
        min_base_quality = config["lcr-modules"]["cfDNA_umi_workflow"].get("filter_min_base_quality", 20),
        max_base_error_rate = config["lcr-modules"]["cfDNA_umi_workflow"].get("filter_max_base_error_rate", 0.3),
        consensus_threads = 4,
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    threads: 2
    resources:
        mem_mb = fgbio_mem
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.call_and_filter_consensus.log")
    shell:
        """
        fgbio -Xmx6g --compression 0 CallDuplexConsensusReads \
            --input {input.bam} \
            --output /dev/stdout \
            --min-reads {params.min_reads} \
            --min-input-base-quality {params.min_input_base_quality} \
            --threads {params.consensus_threads} 2>> {log} \
        | fgbio -Xmx4g --compression 1 FilterConsensusReads \
            --input /dev/stdin \
            --output {output.bam} \
            --ref {params.ref_w_dict} \
            --min-reads {params.min_reads} \
            --min-base-quality {params.min_base_quality} \
            --max-base-error-rate {params.max_base_error_rate} 2>> {log}
        """

# =============================================================================
# STEP 3: Consensus Unmapped BAM -> Final BAM
#         (realign, remove tags, sort, index)
# =============================================================================

rule realign_and_finalize:
    """
    Re-align consensus reads, remove bulky tags, sort, and index.
    
    Pipeline:
      1. samtools fastq: Convert consensus unmapped BAM to FASTQ
      2. bwa mem: Align consensus reads
      3. fgbio ZipperBams: Merge with unmapped BAM tags
         --tags-to-reverse/revcomp Consensus: Keep per-base tags in sync
      4. samtools view -x: Remove bulky per-base tags (ad, ac, bd, bc, aq, bq, ae, be)
      5. samtools sort: Coordinate sort with index
    
    Tags preserved: cD (consensus depth), aD, bD, RX (UMI), MI (molecule ID)
    Tags removed: ad, ac, bd, bc, aq, bq, ae, be (per-base arrays, large)
    """
    input:
        unmapped_bam = rules.call_and_filter_consensus.output.bam,
    output:
        bam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.tempestv2.consensus.bam"),
    params:
        sort_threads = 4,
        # genome references
        bwa_ref = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"],
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["bwa_threads"]
    resources:
        mem_mb = bwa_mem
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.realign_finalize.log")
    shell: 
        """
    samtools fastq {input.unmapped_bam} 2>> {log} \
    | bwa mem -t {threads} -p -K 150000000 -Y {params.bwa_ref} - 2>> {log} \
    | fgbio -Xmx4g --compression 0 --async-io ZipperBams \
        --unmapped {input.unmapped_bam} \
        --ref {params.ref_w_dict} \
        --tags-to-reverse Consensus \
        --tags-to-revcomp Consensus 2>> {log} \
    | samtools view -h \
        -x ad -x ac -x bd -x bc -x aq -x bq -x ae -x be 2>> {log} \
    | samtools sort -@ {params.sort_threads} -o {output.bam} 2>> {log} && \
     samtools index {output.bam} 2>> {log}
        """

# =============================================================================
# QUALITY CONTROL RULES
# =============================================================================

rule qc_fastqc:
    """Run FastQC on unmapped BAM (raw reads with UMIs)."""
    input:
        bam = rules.fgbio_fastq_to_bam.output.bam
    output:
        html = os.path.join(BAM_OUTDIR, "Q1-fastqc", "{sample}.unmapped_fastqc.html")
    params:
        outdir = os.path.join(BAM_OUTDIR, "Q1-fastqc")
    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.fastqc.log")
    shell:
        "fastqc -o {params.outdir} --nogroup -f bam {input.bam} &> {log}"

rule qc_picard_hsmetrics:
    """Collect hybrid selection metrics on final consensus BAM."""
    input:
        bam = rules.realign_and_finalize.output.bam,
    output:
        hsmet = os.path.join(BAM_OUTDIR, "Q2-hs_metrics", "{sample}.hs_metrics.txt"),
        tarcov = os.path.join(BAM_OUTDIR, "Q2-hs_metrics", "{sample}.target_coverage.txt")
    params:
        capture_reg_il = lambda wildcards: get_capture_regions(wildcards),
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    threads: 4
    resources:
        mem_mb = 5000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.hsmetrics.log")
    shell: 
        """
    picard CollectHsMetrics \
        -R {params.ref_w_dict} \
        -TI {params.capture_reg_il} \
        -BI {params.capture_reg_il} \
        -I {input.bam} \
        -O {output.hsmet} \
        --PER_TARGET_COVERAGE {output.tarcov} \
        --MAX_RECORDS_IN_RAM 5000000 \
        --COVERAGE_CAP 20000 &> {log}
        """

rule qc_picard_oxog:
    """Collect OxoG metrics (8-oxoguanine artifacts)."""
    input:
        bam = rules.realign_and_finalize.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR, "Q3-oxog_metrics", "{sample}.oxoG_metrics.txt")
    params:
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]

    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.oxog.log")
    shell:
        "picard CollectOxoGMetrics -I {input.bam} -R {params.ref_w_dict} -O {output.txt} &> {log}"

rule qc_picard_insertsize:
    """Collect insert size metrics."""
    input:
        bam = rules.realign_and_finalize.output.bam
    output:
        txt = os.path.join(BAM_OUTDIR, "Q4-insert_size", "{sample}.insert_size_metrics.txt"),
        pdf = os.path.join(BAM_OUTDIR, "Q4-insert_size", "{sample}.insert_size_histogram.pdf")
    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.insertsize.log")
    shell:
        "picard CollectInsertSizeMetrics -I {input.bam} -O {output.txt} -H {output.pdf} &> {log}"

rule qc_fgbio_errorrate:
    """Calculate error rate by read position."""
    input:
        bam = rules.realign_and_finalize.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR, "Q5-error_rate", "{sample}.error_rate_by_read_position.txt")
    params:
        outprefix = os.path.join(BAM_OUTDIR, "Q5-error_rate", "{sample}"),
        capture_reg_il = lambda wildcards: get_capture_regions(wildcards),
        dbsnpvcf = config["lcr-modules"]["cfDNA_umi_workflow"]["dbsnpvcf"],
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.errorrate.log")
    shell:
        """
    fgbio ErrorRateByReadPosition \
        -i {input.bam} \
        -r {params.ref_w_dict} \
        -v {params.dbsnpvcf} \
        -l {params.capture_reg_il} \
        --collapse \
        -o {params.outprefix} &> {log}
        """

rule qc_calc_dupl:
    """Calculate duplication/collapse rate (raw reads vs consensus)."""
    input:
        collapsed_bam = rules.realign_and_finalize.output.bam,
        all_reads_bam = rules.fgbio_fastq_to_bam.output.bam
    output:
        txt = os.path.join(BAM_OUTDIR, "Q6-dupl_rate", "{sample}.dup_metrics.txt")
    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/pysam.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.calc_dupl.log")
    params:
        script = os.path.join(UTILSDIR, "qc_calc_dupl.py")
    shell: 
        """
    python {params.script} \
    --collapsed_bam {input.collapsed_bam} \
    --all_reads_bam {input.all_reads_bam} \
    --output {output.txt} \
    --sample {wildcards.sample} &> {log}
        """

rule qc_validate_sam:
    """Validate final BAM file."""
    input:
        bam = rules.realign_and_finalize.output.bam,
    output:
        txt = os.path.join(BAM_OUTDIR, "Q7-validatesam", "{sample}.ValidateSamFile.txt")
    params:
        ref_w_dict = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    threads: 1
    resources:
        mem_mb = 5000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.validatesam.log")
    shell:
        "picard ValidateSamFile -I {input.bam} -R {params.ref_w_dict} > {output.txt} 2> {log}"

# =============================================================================
# TARGET RULES
# =============================================================================

rule all_bams:
    """Generate all final BAMs and QC outputs."""
    input:
        # Final BAMs
        expand(rules.realign_and_finalize.output.bam, sample=UMI_SAMPLESHEET["sample_id"]),
        # QC outputs
        expand(rules.qc_fastqc.output.html, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_picard_hsmetrics.output.hsmet, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_picard_oxog.output.txt, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_picard_insertsize.output.txt, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_fgbio_errorrate.output.txt, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_calc_dupl.output.txt, sample=UMI_SAMPLESHEET["sample_id"]),
        expand(rules.qc_validate_sam.output.txt, sample=UMI_SAMPLESHEET["sample_id"])
