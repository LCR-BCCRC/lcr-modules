#!/usr/bin/env snakemake

import os
import sys
import glob
import pandas as pd

sys.path.append(os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "/modules/tempest/1.0/")) # add local module to path
import utils.fastq_utils as fu

# generate paths for file locations
BAM_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "bam_pipeline")
UTILSDIR = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "/modules/tempest/1.0/utils")

UMI_SAMPLESHEET = config["lcr-modules"]["_shared"]["samples"].copy()

# Helper function to get capture regions for a sample
def get_capture_regions(wildcards):
    """Get the capture regions file path for a sample from the sample sheet and config"""
    capture_space = UMI_SAMPLESHEET.loc[UMI_SAMPLESHEET["sample_id"] == wildcards.sample, "capture_space"]
    
    # check if capture space is found in sample sheet
    if capture_space.empty:
        raise ValueError(f"Capture space not found for {wildcards.sample}")
    
    capture_space_value = capture_space.values[0]
    
    # check if capture space is found in config
    if capture_space_value not in config["lcr-modules"]["_shared"]["captureregionsil"]:
        raise ValueError(f"Capture space '{capture_space_value}' not found in config")
        
    return config["lcr-modules"]["_shared"]["captureregionsil"][capture_space_value]

rule trim_umi:
    """
    Trim the UMI from the front (and end) of each read pair
    Also does some FASTQ cleanup (polyG tails etc)

    Note that we are using --trim_front here instead of --umi, as we do not want the UMI stored in the read name
    as then we can't use fgbio AnnotateBamWithUMIs, as the read names will be different from the FASTQ file
    (In theory we could use CopyUmiFromReadName, but that requires a "-" deliminator between the forward and reverse
    UMIs while fastp uses a "_" so)
    """
    input:
        r1 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r1_path"],
        r2 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r2_path"]
    output:
        r1 = temp(os.path.join(BAM_OUTDIR, "01-trimmedfastqs", "{sample}.R1.trimmed.fastq.gz")),
        r2 = temp(os.path.join(BAM_OUTDIR, "01-trimmedfastqs", "{sample}.R2.trimmed.fastq.gz")),
        fastp_report = os.path.join(BAM_OUTDIR, "01-trimmedfastqs", "{sample}.fastp.json")
    params:
        barcodelength = config["lcr-modules"]["cfDNA_umi_workflow"]["barcodelength"],
        outdir = os.path.join(BAM_OUTDIR, "01-trimmedfastqs")
    threads: 4
    resources:
        mem_mb = 10000
    conda:
        "envs/fastp.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs", "{sample}.fastp.log")
    shell: """
fastp --overrepresentation_analysis --detect_adapter_for_pe --trim_front1 {params.barcodelength} \
--trim_front2 {params.barcodelength} --in1 {input.r1} --in2 {input.r2} --thread {threads} --out1 {output.r1} --out2 {output.r2} \
--report_title "fastp {wildcards.sample}" --json {output.fastp_report} --trim_poly_x \
--qualified_quality_phred 20 2> {log}
"""

rule bwa_align_unsorted:
    input:
        r1 = rules.trim_umi.output.r1,
        r2 = rules.trim_umi.output.r2,
        r1bam = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r1_path"],
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        bam = temp(os.path.join(BAM_OUTDIR , "02-BWA" , "{sample}.bwa.unsort.bam"))
    params:
        readgroup = lambda w, input: fu.generate_read_group(input.r1bam, w.sample, config),
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["bwa_threads"]
    resources:
        mem_mb = 10000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs" , "{sample}.bwa_allreads.log")
    shell:
        "bwa mem -t {threads} -R \"{params.readgroup}\" {input.refgenome} {input.r1} {input.r2} 2> {log} | samtools view -b > {output.bam} 2>> {log}"

def ann_umi_mem(wildcards, attempt, input):
    if attempt == 1:
        return max(35000, 6 * input.size_mb)
    elif attempt == 2:
        return (input.size_mb + 35000)
    elif attempt == 3:
        return (input.size_mb + 65000)

# Add UMI tag, mem hungry step
rule fgbio_annotate_umis:
    input:
        bam = rules.bwa_align_unsorted.output.bam,
        r1 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r1_path"],
        r2 = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r2_path"]
    output:
        bam = temp(os.path.join(BAM_OUTDIR , "03-withumis" , "{sample}.bwa.umi.namesort.bam"))
    params:
        umiloc = config["lcr-modules"]["cfDNA_umi_workflow"]["barcodelocation"],
        java_temp = config["lcr-modules"]["cfDNA_umi_workflow"]["java_temp_dir"],
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.annotateumis.log")
    resources:
        mem_mb = ann_umi_mem
    threads: 12
    shell:
        "fgbio -Xms500m -Xmx{resources.mem_mb}m -Djava.io.tmpdir={params.java_temp} AnnotateBamWithUmis --input {input.bam} --fastq {input.r1} --fastq {input.r2} --read-structure {params.umiloc} --output {output.bam} &> {log}"

# Group reads by UMI into families
rule fgbio_group_umis:
    input:
        bam = rules.fgbio_annotate_umis.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        bam = os.path.join(BAM_OUTDIR, "04-umigrouped", "{sample}.umigrouped.sort.bam"),
        txt = os.path.join(BAM_OUTDIR , "04-umigrouped" , "{sample}.umigrouped.famsize.txt")
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["samtools_sort_threads"]
    resources:
        mem_mb = ann_umi_mem
    params:
        maxedits = config["lcr-modules"]["cfDNA_umi_workflow"]["umiedits"],
        outdir = os.path.join(BAM_OUTDIR ,"04-umigrouped"),
        mem_per_t =  lambda wildcards, threads, resources: int(resources.mem_mb // threads),
        java_temp = config["lcr-modules"]["cfDNA_umi_workflow"]["java_temp_dir"],
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR , "logs", "{sample}.groupumis.log")
    shell:
        """samtools sort -m {params.mem_per_t}M -@ {threads} -n {input.bam} | 
        fgbio -Djava.io.tmpdir={params.java_temp} SetMateInformation --ref {input.refgenome} 2> {log} | 
        fgbio -Djava.io.tmpdir={params.java_temp} GroupReadsByUmi --edits {params.maxedits} --family-size-histogram {output.txt} --strategy paired > {output.bam} 2>> {log}"""

# Generate a consensus of these families
rule fgbio_duplex_consensus:
    input:
        bam = rules.fgbio_group_umis.output.bam
    output:
        bam = temp(os.path.join(BAM_OUTDIR , "05-duplexconsensus" , "{sample}.consensus.unmapped.bam"))
    params:
        minreads = config["lcr-modules"]["cfDNA_umi_workflow"]["minreads"],
        sampleid = "{sample}"
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["duplexconsensus_threads"]
    resources:
        mem_mb = 10000
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.duplexconsensus.log")
    shell:
        "fgbio CallDuplexConsensusReads --input {input.bam} --output {output.bam} --threads {threads} --read-group-id {params.sampleid} --min-reads {params.minreads} &> {log}"


# Because CallDuplexConsensusReads recalculates base qualities, and those numbers can be above those supported by certain tools, set a upper
# limit to base qualities.
# Also, remove all the extra space-taking tags
rule sanitize_bam:
    input:
        bam = rules.fgbio_duplex_consensus.output.bam
    output:
        bam = temp(os.path.join(BAM_OUTDIR , "06-sanitizebam" , "{sample}.consensus.unmapped.capqual.bam"))
    params:
        max_base_qual = int(config["lcr-modules"]["cfDNA_umi_workflow"]["max_base_qual"]),  # Bases with quality scores above this are capped at this
        tagstoremove = config["lcr-modules"]["cfDNA_umi_workflow"]["tags_to_remove"],
        min_base_qual = int(config["lcr-modules"]["cfDNA_umi_workflow"]["min_base_qual"])  # Bases with quality scores below this are masked
    conda:
        "envs/pysam.yaml"
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["basequal_threads"]
    resources:
        mem_mb = 5000
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.sanitizebam.log")
    shell:
        f"""python {os.path.join(UTILSDIR, "sanitize_bam.py")} --in_bam {{input.bam}} --out_bam {{output.bam}} \
        --tagstoremove {{params.tagstoremove}} --max_base_qual {{params.max_base_qual}} \
        --threads {{threads}} &> {{log}}"""

# Covert unaligned BAM back to FASTQ and map reads
rule bwa_realign_bam:
    input:
        bam = rules.sanitize_bam.output.bam,
        r1bam = config["lcr-modules"]["cfDNA_umi_workflow"]["fastq_r1_path"],
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        bam = temp(os.path.join(BAM_OUTDIR, "07-consensus_aligned", "{sample}.consensus.mapped.namesort.bam"))
    threads:
        config["lcr-modules"]["cfDNA_umi_workflow"]["bwa_threads"]
    resources:
        mem_mb = 15000
    params:
        readgroup = lambda w, input: fu.generate_read_group(input.r1bam, w.sample, config)
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR , "logs", "{sample}.bwa_realign.log")
    shell:
        """samtools fastq -@ {threads} -N {input.bam} 2> {log} | bwa mem -R \"{params.readgroup}\" -p -t {threads} {input.refgenome} - 2>> {log} | samtools view -b |
        picard SortSam -SO queryname -I /dev/stdin -O /dev/stdout 2>> {log} | fgbio SetMateInformation --ref {input.refgenome} --output {output.bam} &> {log}"""

def pic_ann_mem(wildcards, input,attempt):
    if attempt == 1:
        return max(8 * input.size_mb, 10000)
    elif attempt == 2:
        return max(15 * input.size_mb, 10000)
    elif attempt == 3:
        return max(30 * input.size_mb, 10000)
    else:
        return max(50 * input.size_mb, 10000)

# Add back in family information from the unaligned consensus BAM
rule picard_annotate_bam:
    input:
        unaligned_bam = rules.sanitize_bam.output.bam,
        aligned_bam = rules.bwa_realign_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    output:
        bam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.picardannotate.log")
    resources:
        # make adapttive mem requirements
        # mem_mb= lambda wc,input, attempt: pic_ann_mem(wc, input,attempt)
        mem_mb = pic_ann_mem
    threads: 4
    shell:
        """picard -Xms500m -Xmx{resources.mem_mb}m SortSam -I {input.unaligned_bam} -O /dev/stdout -SO queryname --REFERENCE_SEQUENCE {input.refgenome} |
           picard -Xms500m -Xmx10g MergeBamAlignment --ALIGNED_BAM {input.aligned_bam} --UNMAPPED_BAM /dev/stdin --REFERENCE_SEQUENCE {input.refgenome} -O {output.bam} &> {log} && 
           samtools index {output.bam}"""

# Output sentinel confirming that the final BAMs are valid
rule picard_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR , "09-validoutput" , "{sample}.consensus.mapped.ValidateSamFile.is_valid")
    params:
        outdir = os.path.join(BAM_OUTDIR , "09-validoutput")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR, "logs" , "{sample}.picardvalidatesam.log")
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"


### QUALITY CONTROL RULES ###
### Currently supported tools:
### FASTQC
### PICARD

rule qc_fastqc:
    input:
        bam = rules.bwa_align_unsorted.output.bam
    output:
        qc = os.path.join(BAM_OUTDIR , "Q1-fastqc" , "{sample}.bwa.unsort_fastqc.html")
    params:
        outdir = os.path.join(BAM_OUTDIR ,"Q1-fastqc")
    conda:
        "envs/picard_fastqc.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.fastqc.log")
    shell:
        "fastqc -o {params.outdir} --nogroup -f bam {input.bam} 2> {log}"

rule qc_picard_hsmetrics:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        hsmet = os.path.join(BAM_OUTDIR , "Q2-hs_metrics" , "{sample}.hs_metrics.txt"),
        tarcov = os.path.join(BAM_OUTDIR , "Q2-hs_metrics" , "{sample}.target_coverage.txt")
    params:
        capture_reg_il = lambda wildcards: get_capture_regions(wildcards),
        outdir = os.path.join(BAM_OUTDIR , "Q2-hs_metrics"),
        max_ram_records = "5000000",
        cov_cap_sens = "20000"
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 4
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.picard_hsmet.log")
    shell:
        "picard CollectHsMetrics -R {input.refgenome} -TI {params.capture_reg_il} -BI {params.capture_reg_il} -I {input.bam} -O {output.hsmet} --PER_TARGET_COVERAGE {output.tarcov} --MAX_RECORDS_IN_RAM {params.max_ram_records} --COVERAGE_CAP {params.cov_cap_sens} 2> {log}"

rule qc_picard_oxog:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR, "Q3-oxog_metrics", "{sample}.oxoG_metrics.txt")
    params:
        outdir = os.path.join(BAM_OUTDIR , "Q3-oxog_metrics")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR , "logs", "{sample}.picard_oxoG.log")
    shell:
        "picard CollectOxoGMetrics -I {input.bam} -R {input.refgenome} -O {output.txt} 2> {log}"

rule qc_picard_insertsize:
    input:
        bam = rules.picard_annotate_bam.output.bam
    output:
        txt = os.path.join(BAM_OUTDIR, "Q4-insert_size", "{sample}.insert_size_metrics.txt"),
        pdf = os.path.join(BAM_OUTDIR , "Q4-insert_size" , "{sample}.insert_size_histogram.pdf")
    params:
        outdir = os.path.join(BAM_OUTDIR , "Q4-insert_size")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.picard_insertsize.log")
    shell:
        "picard CollectInsertSizeMetrics -I {input.bam} -O {output.txt} -H {output.pdf} 2> {log}"

rule qc_fgbio_errorrate:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["picard_refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR , "Q5-error_rate" , "{sample}.error_rate_by_read_position.txt")
    params:
        outprefix = os.path.join(BAM_OUTDIR, "Q5-error_rate" , "{sample}"),
        outdir = os.path.join(BAM_OUTDIR , "Q5-error_rate"),
        capture_reg_il = lambda wildcards: get_capture_regions(wildcards),
        dbsnpvcf = config["lcr-modules"]["cfDNA_umi_workflow"]["dbsnpvcf"]
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.error_rate_by_position.log")
    shell:
        "fgbio ErrorRateByReadPosition -i {input.bam} -r {input.refgenome} -v {params.dbsnpvcf} -l {params.capture_reg_il} --collapse -o {params.outprefix} 2> {log}"

rule qc_calc_dupl:
    input:
        collapsed_bam = rules.picard_annotate_bam.output.bam,
        all_reads_bam = rules.bwa_align_unsorted.output.bam
    output:
        txt = os.path.join(BAM_OUTDIR , "Q6-dupl_rate" , "{sample}.dup_metrics.txt")
    params:
        outdir = os.path.join(BAM_OUTDIR, "Q6-dupl_rate")
    conda:
        "envs/pysam.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        f"""python {os.path.join(UTILSDIR, "qc_calc_dupl.py")} --collapsed_bam {{input.collapsed_bam}} --all_reads_bam {{input.all_reads_bam}} \
        --output {{output.txt}} --sample {{wildcards.sample}}"""

# Output sentinel confirming that the final BAMs are valid
rule qc_validate_sam:
    input:
        bam = rules.picard_annotate_bam.output.bam,
        refgenome = config["lcr-modules"]["cfDNA_umi_workflow"]["refgenome"]
    output:
        txt = os.path.join(BAM_OUTDIR , "Q7-validatesam" , "{sample}.consensus.mapped.ValidateSamFile.is_valid")
    params:
        outdir = os.path.join(BAM_OUTDIR , "Q7-validatesam")
    conda:
        "envs/bwa_picard_fgbio.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(BAM_OUTDIR , "logs" , "{sample}.picardvalidatesam.log")
    shell:
        "picard ValidateSamFile -I {input.bam} -R {input.refgenome} > {output.txt} 2> {log}"

# Merge QC results via multiqc
checkpoint qc_multiqc:
    input:
        # Run multiqc once per run, and merge all samples from that run
        dupl = lambda w: list(os.path.join(BAM_OUTDIR, "Q6-dupl_rate", x + ".dup_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        errorate = lambda w: list(os.path.join(BAM_OUTDIR, "Q5-error_rate", x + ".error_rate_by_read_position.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        insertsize = lambda w: list(os.path.join(BAM_OUTDIR,"Q4-insert_size", x + ".insert_size_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        oxog = lambda w: list(os.path.join(BAM_OUTDIR, "Q3-oxog_metrics", x + ".oxoG_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        hsmet = lambda w: list(os.path.join(BAM_OUTDIR, "Q2-hs_metrics", x + ".hs_metrics.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        fastp = lambda w: list(os.path.join(BAM_OUTDIR, "01-trimmedfastqs", x + ".fastp.json") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        fastqc = lambda w: list(os.path.join(BAM_OUTDIR, "Q1-fastqc", x + ".bwa.unsort_fastqc.html") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        validatesam = lambda w: list(os.path.join(BAM_OUTDIR, "Q7-validatesam", x + ".consensus.mapped.ValidateSamFile.is_valid") for x in SAMPLELIST if sample_to_runid[x] == w.runid),
        famsizehist = lambda w: list(os.path.join(BAM_OUTDIR, "04-umigrouped", x + ".umigrouped.famsize.txt") for x in SAMPLELIST if sample_to_runid[x] == w.runid)
    output:
        html = os.path.join(BAM_OUTDIR , "Q9-multiqc" , "multiqc_report.{runid}.html"),
    params:
        outdir = os.path.join(BAM_OUTDIR , "Q9-multiqc"),
        outname = lambda w: "multiqc_report." + w.runid + ".html",
        ignoresamples = lambda w: "\' --ignore-samples \'".join(x for x in SAMPLELIST if sample_to_runid[x] != w.runid),
        modules = "-m picard -m fastqc -m fgbio -m fastp",  # Should start with -m flag
        config = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "config/multiqc_config.yaml"),
        dupl_dir = rules.qc_calc_dupl.params.outdir,
        errorrate_dir = rules.qc_fgbio_errorrate.params.outdir,
        insertsize_dir = rules.qc_picard_insertsize.params.outdir,
        oxog_dir = rules.qc_picard_oxog.params.outdir,
        hsmet_dir = rules.qc_picard_hsmetrics.params.outdir,
        fastp_dir = rules.trim_umi.params.outdir,
        fastqc_dir = rules.qc_fastqc.params.outdir,
        validsam_dir = rules.qc_validate_sam.params.outdir,
        famsize_dir = rules.fgbio_group_umis.params.outdir
    conda:
        "envs/picard_fastqc.yaml"
    log:
        os.path.join(BAM_OUTDIR, "logs" , "multiqc_{runid}.log")
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        """multiqc --no-data-dir --interactive --config {params.config} --outdir {params.outdir} --filename {params.outname} \
        --force {params.modules} {params.dupl_dir} {params.errorrate_dir} {params.insertsize_dir} {params.oxog_dir} {params.hsmet_dir} {params.fastqc_dir} {params.validsam_dir} {params.famsize_dir} {params.fastp_dir} \
        --ignore-samples \'{params.ignoresamples}\' > {log}"""


rule all_bams:
    input:
        expand(str(rules.qc_calc_dupl.output.txt), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_fgbio_errorrate.output.txt), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_picard_insertsize.output.txt), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_picard_oxog.output.txt), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_picard_hsmetrics.output.hsmet), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_fastqc.output.qc), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.qc_validate_sam.output.txt), sample=UMI_SAMPLESHEET["sample_id"]),
        expand(str(rules.fgbio_group_umis.output.txt), sample=UMI_SAMPLESHEET["sample_id"])
