#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Brett Collinge
# Module Author:    Brett Collinge
# Contributors:     N/A


##### SETUP #####


import oncopipe as op
import pandas as pd
import os

# Check oncopipe version (consistent with other modules in this project)
min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
        '\x1b[0;31;40m' + f"ERROR: oncopipe {current_version} < required {min_oncopipe_version}" + '\x1b[0m'
    )
    sys.exit("Please update oncopipe in your environment")


# Setup module and store module-specific configuration in `CFG`. This now lives in
# lcr-modules/modules/gatk_cnv/1.0, so op.setup_module resolves it normally (the
# cs-modules path-swap hack from the original crispr_screens version is removed).
CFG = op.setup_module(
    name = "gatk_cnv",
    version = "1.0",
    subdirectories = ["inputs", "intervals", "pon", "counts", "denoised", "model", "called", "purecn", "gene_matrix", "mutect2", "outputs"]
)


# Capture data that op.cleanup_module will pop from CFG (samples, runs, paired_runs,
# unpaired_runs) into module globals so deferred-eval input functions can still
# reach them after cleanup. The other CFG sub-dicts (dirs, inputs, options) are
# not popped — they remain reachable via config["lcr-modules"]["gatk_cnv"][...].
_RUNS = CFG["runs"]


# Load the PoN sample table directly from the TSV path. We deliberately don't
# cache the DataFrame in CFG — op.cleanup_module yaml-dumps the whole CFG and
# fails on non-serializable objects nested below the well-known top-level keys.
_pon_tsv = CFG["pon"].get("metadata_tsv")
if not _pon_tsv or _pon_tsv == "__UPDATE__":
    sys.exit("gatk_cnv: pon.metadata_tsv must be set in config/gatk_cnv.yaml")
PON_SAMPLES = pd.read_csv(_pon_tsv, sep="\t", dtype=str)
_pon_path_col = CFG["pon"]["bam_path_column"]

if "sample_id" not in PON_SAMPLES.columns or _pon_path_col not in PON_SAMPLES.columns:
    sys.exit(f"gatk_cnv: PoN metadata must have 'sample_id' and '{_pon_path_col}' columns.")

if len(PON_SAMPLES) < CFG["pon"]["min_samples"]:
    sys.exit(
        f"gatk_cnv: PoN has {len(PON_SAMPLES)} samples; need at least "
        f"{CFG['pon']['min_samples']}."
    )

# Build a lookup: sample_id -> absolute BAM path (from the PoN metadata)
_PON_BAM_BY_ID = dict(zip(PON_SAMPLES["sample_id"], PON_SAMPLES[_pon_path_col]))


# Capture-space identifiers (resolved against config/reference_files.yaml)
_CASE_CAPSPACE = CFG["capture_spaces"]["case"]
_PON_CAPSPACE  = CFG["capture_spaces"]["pon"]

# Genome build for PoN runs. PoN samples are not in CFG["runs"]; we hardcode grch37
# since that's the only build currently supported by the registered capture_spaces.
_PON_GENOME_BUILD = "grch37"
_PON_SEQ_TYPE     = "capture"

# Force pure-Java (F2J) BLAS for GATK Spark-based tools (CreateReadCountPanelOfNormals,
# DenoiseReadCounts). Cluster compute nodes have an incomplete system BLAS missing
# `cblas_dspr`, so netlib-native fails with `symbol lookup error` when Spark tries to
# do SVD. F2J is the documented netlib fallback — pure Java, ~3-10x slower than native
# but irrelevant for our small matrices (single-digit-eigensample SVD on ~200K intervals).
_F2J_BLAS = (
    "-Dcom.github.fommil.netlib.BLAS=com.github.fommil.netlib.F2jBLAS "
    "-Dcom.github.fommil.netlib.LAPACK=com.github.fommil.netlib.F2jLAPACK "
    "-Dcom.github.fommil.netlib.ARPACK=com.github.fommil.netlib.F2jARPACK"
)


# Define rules to be run locally (no cluster submission)
localrules:
    _gatk_cnv_input_bam,
    _gatk_cnv_input_pon_bam,
    _gatk_cnv_output_called_seg,
    _gatk_cnv_output_denoised,
    _gatk_cnv_output_plots,
    _gatk_cnv_output_purecn,
    _gatk_cnv_output_gene_matrix,
    _gatk_cnv_output_summary,
    _gatk_cnv_all,


##### INPUT SYMLINKS #####


# Symlink case sample BAMs/CRAMs into 00-inputs/. The CRAM .crai index is
# also exposed as both .crai and .bai because some downstream tools look for
# .bai even when given a CRAM.
rule _gatk_cnv_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Symlink PoN sample BAMs/CRAMs into 00-inputs/pon/. PoN bam paths come from the
# metadata TSV; the dz wrapper points them at indexed, stably-named files (CRAM
# stays CRAM). Expose the index as both .bam.bai and .bam.crai so GATK finds it
# whether the underlying file is BAM or CRAM.
rule _gatk_cnv_input_pon_bam:
    input:
        bam = lambda w: _PON_BAM_BY_ID[w.sample_id],
        bai = lambda w: _PON_BAM_BY_ID[w.sample_id] + ".bai"
    output:
        bam = CFG["dirs"]["inputs"] + "pon/bam/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "pon/bam/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "pon/bam/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


##### REFERENCE / INTERVAL PREP #####


# V5 ∩ V6+UTR intersection BED, autosomes only. Both BEDs come from the
# reference_files subworkflow (registered in config/reference_files.yaml).
rule _gatk_cnv_intersect_bed:
    input:
        v5_bed = reference_files("genomes/{genome_build}/capture_space/" + _CASE_CAPSPACE + ".bed"),
        v6_bed = reference_files("genomes/{genome_build}/capture_space/" + _PON_CAPSPACE  + ".bed")
    output:
        bed = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.autosomes.bed"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/intersect_bed.stderr.log"
    params:
        # b37/grch37 BEDs use no "chr" prefix; sex chromosomes are literal "X" and "Y".
        # If include_sex_chromosomes=True, retain X and Y; otherwise drop them.
        keep_sex = CFG["options"]["include_sex_chromosomes"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads: 1
    shell:
        op.as_one_line("""
        bedtools intersect -a {input.v5_bed} -b {input.v6_bed} 2> {log.stderr}
            |
        awk -v keep_sex={params.keep_sex} 'BEGIN {{FS=OFS="\\t"}}
            {{
              chrom = $1;
              gsub(/^chr/, "", chrom);
              if (keep_sex == "True" || (chrom != "X" && chrom != "Y" && chrom != "MT" && chrom != "M")) print $0;
            }}' > {output.bed} 2>> {log.stderr}
        """)


rule _gatk_cnv_preprocess_intervals:
    input:
        bed = str(rules._gatk_cnv_intersect_bed.output.bed),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        ref_dict  = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.preprocessed.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/preprocess_intervals.stderr.log"
    params:
        bin_length = CFG["options"]["bin_length"],
        padding    = CFG["options"]["padding"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["preprocess_intervals"]
    resources:
        **CFG["resources"]["preprocess_intervals"]
    shell:
        op.as_one_line("""
        gatk PreprocessIntervals
            -R {input.ref_fasta}
            -L {input.bed}
            --bin-length {params.bin_length}
            --padding {params.padding}
            --interval-merging-rule OVERLAPPING_ONLY
            -O {output.interval_list}
            2> {log.stderr}
        """)


rule _gatk_cnv_annotate_intervals:
    input:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.case_filtered.preprocessed.interval_list",
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        ref_dict  = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        annotated = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.case_filtered.annotated.tsv"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/annotate_intervals.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["annotate_intervals"]
    resources:
        **CFG["resources"]["annotate_intervals"]
    shell:
        op.as_one_line("""
        gatk AnnotateIntervals
            -R {input.ref_fasta}
            -L {input.interval_list}
            --interval-merging-rule OVERLAPPING_ONLY
            -O {output.annotated}
            2> {log.stderr}
        """)


# Build the common biallelic SNP VCF for CollectAllelicCounts. Restrict gnomAD
# to autosomal intersection intervals + biallelic SNPs at AF>0.001.
rule _gatk_cnv_biallelic_snps:
    input:
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        bed    = str(rules._gatk_cnv_intersect_bed.output.bed)
    output:
        vcf = CFG["dirs"]["intervals"] + "{genome_build}/common_biallelic.vcf.gz",
        tbi = CFG["dirs"]["intervals"] + "{genome_build}/common_biallelic.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/biallelic_snps.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["biallelic_snps"]
    resources:
        **CFG["resources"]["biallelic_snps"]
    shell:
        op.as_one_line("""
        bcftools view -R {input.bed} -m2 -M2 -v snps -i 'INFO/AF>0.001'
            --threads {threads} -Oz -o {output.vcf} {input.gnomad} 2> {log.stderr}
            &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


##### READ COUNTS — CASE SAMPLES #####


# Initial-pass case counts at the UNFILTERED intervals. Used as input to the
# case-coverage filter, which then drives the main pass. This produces ~228k
# intervals' worth of counts; the main collect_counts rule (below) writes the
# ~213k-interval counts that downstream rules actually consume.
rule _gatk_cnv_collect_counts_initial:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        interval_list = str(rules._gatk_cnv_preprocess_intervals.output.interval_list),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        counts = CFG["dirs"]["counts"] + "init/{seq_type}--{genome_build}/{sample_id}.counts.hdf5"
    log:
        stderr = CFG["logs"]["counts"] + "init/{seq_type}--{genome_build}/{sample_id}.collect_counts_initial.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["collect_read_counts"]
    resources:
        **CFG["resources"]["collect_read_counts"]
    shell:
        op.as_one_line("""
        gatk CollectReadCounts
            -I {input.bam}
            -L {input.interval_list}
            -R {input.ref_fasta}
            --interval-merging-rule OVERLAPPING_ONLY
            --format HDF5
            -O {output.counts}
            2> {log.stderr}
        """)


# Aggregate per-interval medians across all case samples and drop intervals
# whose median is below threshold. Targets V5 baits that don't capture cleanly
# across the cohort. Output replaces the unfiltered preprocessed.interval_list
# everywhere downstream (annotate, collect_counts, PoN).
def _gatk_cnv_filter_case_intervals_input(wildcards):
    cfg = config["lcr-modules"]["gatk_cnv"]
    runs = _RUNS[_RUNS["tumour_genome_build"] == wildcards.genome_build]
    pairs = runs[["tumour_sample_id", "tumour_seq_type"]].drop_duplicates()
    counts_files = [
        cfg["dirs"]["counts"] + f"init/{row.tumour_seq_type}--{wildcards.genome_build}/{row.tumour_sample_id}.counts.hdf5"
        for row in pairs.itertuples(index=False)
    ]
    return {
        "interval_list": cfg["dirs"]["intervals"] + f"{wildcards.genome_build}/v5_intersect_v6.preprocessed.interval_list",
        "counts": counts_files,
        "filter_script": cfg["inputs"]["filter_intervals"],
    }


rule _gatk_cnv_filter_case_intervals:
    input:
        unpack(_gatk_cnv_filter_case_intervals_input)
    output:
        interval_list = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.case_filtered.preprocessed.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/filter_case_intervals.stderr.log"
    params:
        min_median = CFG["options"]["case_coverage_filter_min_median"]
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["filter_case_intervals"]
    resources:
        **CFG["resources"]["filter_case_intervals"]
    shell:
        op.as_one_line("""
        python {input.filter_script}
            --interval_list {input.interval_list}
            --counts {input.counts}
            --output {output.interval_list}
            --min_median {params.min_median}
            2> {log.stderr}
        """)


rule _gatk_cnv_collect_counts:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        interval_list = str(rules._gatk_cnv_filter_case_intervals.output.interval_list),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        counts = CFG["dirs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.counts.hdf5"
    log:
        stderr = CFG["logs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.collect_counts.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["collect_read_counts"]
    resources:
        **CFG["resources"]["collect_read_counts"]
    shell:
        op.as_one_line("""
        gatk CollectReadCounts
            -I {input.bam}
            -L {input.interval_list}
            -R {input.ref_fasta}
            --interval-merging-rule OVERLAPPING_ONLY
            --format HDF5
            -O {output.counts}
            2> {log.stderr}
        """)


rule _gatk_cnv_collect_allelic_counts:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        snps = str(rules._gatk_cnv_biallelic_snps.output.vcf),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        allelic = CFG["dirs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.allelicCounts.tsv"
    log:
        stderr = CFG["logs"]["counts"] + "{seq_type}--{genome_build}/{sample_id}.collect_allelic_counts.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["collect_allelic_counts"]
    resources:
        **CFG["resources"]["collect_allelic_counts"]
    shell:
        op.as_one_line("""
        gatk CollectAllelicCounts
            -I {input.bam}
            -L {input.snps}
            -R {input.ref_fasta}
            -O {output.allelic}
            2> {log.stderr}
        """)


##### READ COUNTS — PoN SAMPLES #####


rule _gatk_cnv_collect_counts_pon:
    input:
        bam = CFG["dirs"]["inputs"] + "pon/bam/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "pon/bam/{sample_id}.bam.bai",
        interval_list = CFG["dirs"]["intervals"] + _PON_GENOME_BUILD + "/v5_intersect_v6.case_filtered.preprocessed.interval_list",
        ref_fasta = reference_files("genomes/" + _PON_GENOME_BUILD + "/genome_fasta/genome.fa")
    output:
        counts = CFG["dirs"]["counts"] + "pon/" + _PON_GENOME_BUILD + "/{sample_id}.counts.hdf5"
    log:
        stderr = CFG["logs"]["counts"] + "pon/" + _PON_GENOME_BUILD + "/{sample_id}.collect_counts.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["collect_read_counts"]
    resources:
        **CFG["resources"]["collect_read_counts"]
    shell:
        op.as_one_line("""
        gatk CollectReadCounts
            -I {input.bam}
            -L {input.interval_list}
            -R {input.ref_fasta}
            --interval-merging-rule OVERLAPPING_ONLY
            --format HDF5
            -O {output.counts}
            2> {log.stderr}
        """)


rule _gatk_cnv_create_pon:
    input:
        counts = expand(
            CFG["dirs"]["counts"] + "pon/" + _PON_GENOME_BUILD + "/{sample_id}.counts.hdf5",
            sample_id = list(PON_SAMPLES["sample_id"])
        ),
        annotated = CFG["dirs"]["intervals"] + _PON_GENOME_BUILD + "/v5_intersect_v6.case_filtered.annotated.tsv"
    output:
        pon = CFG["dirs"]["pon"] + _PON_GENOME_BUILD + "/v5_intersect_v6.cnv.pon.hdf5"
    log:
        stderr = CFG["logs"]["pon"] + _PON_GENOME_BUILD + "/create_pon.stderr.log"
    params:
        counts_args = lambda w, input: " ".join(f"-I {p}" for p in input.counts),
        # Leave ~20% headroom for JVM metaspace, native code, thread stacks.
        # -Xmx == total mem_mb causes JVM to fail allocating metaspace at startup.
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["create_pon"]
    resources:
        **CFG["resources"]["create_pon"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" CreateReadCountPanelOfNormals
            {params.counts_args}
            --annotated-intervals {input.annotated}
            -O {output.pon}
            2> {log.stderr}
        """ % _F2J_BLAS)


##### CASE-SAMPLE CALLING (paired or unpaired) #####


rule _gatk_cnv_denoise_read_counts:
    input:
        counts = CFG["dirs"]["counts"] + "{seq_type}--{genome_build}/{tumour_id}.counts.hdf5",
        pon    = CFG["dirs"]["pon"] + _PON_GENOME_BUILD + "/v5_intersect_v6.cnv.pon.hdf5"
    output:
        std_cr      = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.standardizedCR.tsv",
        denoised_cr = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoisedCR.tsv"
    log:
        stderr = CFG["logs"]["denoised"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoise.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["denoise_read_counts"]
    resources:
        **CFG["resources"]["denoise_read_counts"]
    shell:
        op.as_one_line("""
        gatk --java-options "%s" DenoiseReadCounts
            -I {input.counts}
            --count-panel-of-normals {input.pon}
            --standardized-copy-ratios {output.std_cr}
            --denoised-copy-ratios    {output.denoised_cr}
            2> {log.stderr}
        """ % _F2J_BLAS)


def _gatk_cnv_model_segments_inputs(wildcards):
    """Build inputs for ModelSegments. Always supplies tumour allelic counts;
    supplies normal allelic counts only when pair_status == 'matched'.
    Re-fetches CFG from config because op.cleanup_module deletes the local CFG."""
    cfg = config["lcr-modules"]["gatk_cnv"]
    inp = {
        "denoised_cr": cfg["dirs"]["denoised"] + f"{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.tumour_id}--{wildcards.normal_id}--{wildcards.pair_status}.denoisedCR.tsv",
        "tumour_allelic": cfg["dirs"]["counts"] + f"{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.tumour_id}.allelicCounts.tsv",
    }
    if wildcards.pair_status == "matched":
        inp["normal_allelic"] = cfg["dirs"]["counts"] + f"{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.normal_id}.allelicCounts.tsv"
    return inp


rule _gatk_cnv_model_segments:
    input:
        unpack(_gatk_cnv_model_segments_inputs)
    output:
        seg_modelFinal = CFG["dirs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.seg",
        seg_cr         = CFG["dirs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cr.seg",
        hets           = CFG["dirs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.hets.tsv"
    log:
        stderr = CFG["logs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.model.stderr.log"
    params:
        out_dir = CFG["dirs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}",
        normal_arg = lambda w, input: (
            f"--normal-allelic-counts {input.normal_allelic}" if w.pair_status == "matched" else ""
        ),
        smoothing_cr = CFG["options"]["model_segments_smoothing_cr"],
        smoothing_af = CFG["options"]["model_segments_smoothing_af"],
        # Leave ~20% headroom for JVM metaspace + thread stacks; -Xmx == total
        # mem_mb causes the JVM to fail allocating metaspace at startup.
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["model_segments"]
    resources:
        **CFG["resources"]["model_segments"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" ModelSegments
            --denoised-copy-ratios {input.denoised_cr}
            --allelic-counts {input.tumour_allelic}
            {params.normal_arg}
            --smoothing-credible-interval-threshold-copy-ratio {params.smoothing_cr}
            --smoothing-credible-interval-threshold-allele-fraction {params.smoothing_af}
            --output {params.out_dir}
            --output-prefix {wildcards.tumour_id}
            2> {log.stderr}
        """ % _F2J_BLAS)


rule _gatk_cnv_call_segments:
    input:
        seg_cr = str(rules._gatk_cnv_model_segments.output.seg_cr)
    output:
        called = CFG["dirs"]["called"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.called.seg"
    log:
        stderr = CFG["logs"]["called"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.call.stderr.log"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["call_segments"]
    resources:
        **CFG["resources"]["call_segments"]
    shell:
        op.as_one_line("""
        gatk CallCopyRatioSegments
            --input {input.seg_cr}
            --output {output.called}
            2> {log.stderr}
        """)


rule _gatk_cnv_plot_denoised:
    input:
        std_cr      = str(rules._gatk_cnv_denoise_read_counts.output.std_cr),
        denoised_cr = str(rules._gatk_cnv_denoise_read_counts.output.denoised_cr),
        ref_dict    = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        denoised_plot = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoised.png"
    log:
        stderr = CFG["logs"]["denoised"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.plot_denoised.stderr.log"
    params:
        out_dir = CFG["dirs"]["denoised"] + "{seq_type}--{genome_build}",
        prefix  = "{tumour_id}--{normal_id}--{pair_status}"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["plot_denoised"]
    resources:
        **CFG["resources"]["plot_denoised"]
    shell:
        op.as_one_line("""
        gatk --java-options "%s" PlotDenoisedCopyRatios
            --standardized-copy-ratios {input.std_cr}
            --denoised-copy-ratios     {input.denoised_cr}
            --sequence-dictionary      {input.ref_dict}
            --output                   {params.out_dir}
            --output-prefix            {params.prefix}
            2> {log.stderr}
        """ % _F2J_BLAS)


rule _gatk_cnv_plot_modeled:
    input:
        denoised_cr    = str(rules._gatk_cnv_denoise_read_counts.output.denoised_cr),
        seg_modelFinal = str(rules._gatk_cnv_model_segments.output.seg_modelFinal),
        hets           = str(rules._gatk_cnv_model_segments.output.hets),
        ref_dict       = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        modeled_plot = CFG["dirs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.modeled.png"
    log:
        stderr = CFG["logs"]["model"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.plot_modeled.stderr.log"
    params:
        out_dir = CFG["dirs"]["model"] + "{seq_type}--{genome_build}",
        prefix  = "{tumour_id}--{normal_id}--{pair_status}"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["plot_modeled"]
    resources:
        **CFG["resources"]["plot_modeled"]
    shell:
        op.as_one_line("""
        gatk --java-options "%s" PlotModeledSegments
            --denoised-copy-ratios {input.denoised_cr}
            --allelic-counts       {input.hets}
            --segments             {input.seg_modelFinal}
            --sequence-dictionary  {input.ref_dict}
            --output               {params.out_dir}
            --output-prefix        {params.prefix}
            2> {log.stderr}
        """ % _F2J_BLAS)


##### PureCN INPUTS — Mutect2 (germline+somatic) and IntervalFile #####


# Generate a PureCN-format interval file via PureCN's IntervalFile.R. Output
# has the columns runAbsoluteCN expects: Target, gc_bias, mappability, Gene.
# The Gene column is what makes PureCN's callAlterations() produce per-gene
# absolute CN — the previous synthesized interval file lacked it, so PureCN's
# per-gene calls were always empty.
#
# DepMap uses pre-generated interval files per bait set (e.g. agilent_hg38_intervals.txt).
# We generate ours on-demand from the V5∩V6+UTR intersection BED.
rule _gatk_cnv_purecn_intervals:
    input:
        bed = str(rules._gatk_cnv_intersect_bed.output.bed),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        intervals = CFG["dirs"]["intervals"] + "{genome_build}/v5_intersect_v6.purecn_intervals.txt"
    log:
        stderr = CFG["logs"]["intervals"] + "{genome_build}/purecn_intervals.stderr.log"
    params:
        purecn_genome = lambda w: "hg19" if w.genome_build == "grch37" else ("hg38" if w.genome_build == "grch38" else w.genome_build)
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["purecn_intervals"]
    resources:
        **CFG["resources"]["purecn_intervals"]
    shell:
        op.as_one_line("""
        INTERVAL_SCRIPT=$(Rscript -e
            "cat(system.file('extdata', 'IntervalFile.R', package='PureCN'))"
        )
        &&
        Rscript "$INTERVAL_SCRIPT"
            --in-file {input.bed}
            --fasta {input.ref_fasta}
            --out-file {output.intervals}
            --genome {params.purecn_genome}
            2> {log.stderr}
        """)


# Tumour-only Mutect2 over the full case-filtered intervals with germline-sites
# genotyping. Captures germline het SNPs (BAF~0.5) that PureCN needs to anchor
# segmentation and refine purity. Matches DepMap's invocation:
#   m2_extra_args="--genotype-germline-sites true --genotype-pon-sites true"
# This is independent of slms_3's restricted-positions Mutect2 (which only
# genotypes positions called by Strelka+LoFreq and is somatic-only by design).
rule _gatk_cnv_mutect2_germline:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        ref_dict  = reference_files("genomes/{genome_build}/genome_fasta/genome.dict"),
        gnomad    = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"),
        intervals = str(rules._gatk_cnv_filter_case_intervals.output.interval_list),
        pon       = reference_files("genomes/{genome_build}/gatk/mutect2_pon.{genome_build}.vcf.gz")
    output:
        vcf   = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2.vcf.gz",
        tbi   = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2.vcf.gz.tbi",
        stats = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2.vcf.gz.stats"
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2_germline.stderr.log"
    params:
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["mutect2_germline"]
    resources:
        **CFG["resources"]["mutect2_germline"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" Mutect2
            -R {input.ref_fasta}
            -I {input.bam}
            -L {input.intervals}
            --germline-resource {input.gnomad}
            --panel-of-normals {input.pon}
            --genotype-germline-sites true
            --genotype-pon-sites true
            -O {output.vcf}
            2> {log.stderr}
        """ % _F2J_BLAS)


# FilterMutectCalls labels variants with PASS / germline / etc. in the FILTER
# column. PureCN's filterVcfMuTect uses this column to distinguish somatic
# from germline candidates. We skip the optional read-orientation /
# contamination filters (which target somatic precision); for purity/ploidy
# fitting via PureCN, the basic FilterMutectCalls output is sufficient.
rule _gatk_cnv_filter_mutect_calls:
    input:
        vcf       = str(rules._gatk_cnv_mutect2_germline.output.vcf),
        tbi       = str(rules._gatk_cnv_mutect2_germline.output.tbi),
        stats     = str(rules._gatk_cnv_mutect2_germline.output.stats),
        ref_fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        ref_dict  = reference_files("genomes/{genome_build}/genome_fasta/genome.dict")
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2.filtered.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.mutect2.filtered.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{sample_id}.filter_mutect_calls.stderr.log"
    params:
        java_mem = lambda w, resources: int(resources.mem_mb * 0.8)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["filter_mutect_calls"]
    resources:
        **CFG["resources"]["filter_mutect_calls"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}m %s" FilterMutectCalls
            -R {input.ref_fasta}
            -V {input.vcf}
            --stats {input.stats}
            -O {output.vcf}
            2> {log.stderr}
        """ % _F2J_BLAS)


##### PureCN #####


rule _gatk_cnv_purecn:
    input:
        run_purecn      = ancient(CFG["inputs"]["run_purecn"]),
        seg_modelFinal  = str(rules._gatk_cnv_model_segments.output.seg_modelFinal),
        mutect_vcf      = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{tumour_id}.mutect2.filtered.vcf.gz",
        intervals       = str(rules._gatk_cnv_purecn_intervals.output.intervals),
        ref_fasta       = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        purity_csv = CFG["dirs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.csv",
        genes_csv  = CFG["dirs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_genes.csv",
        loh_csv    = CFG["dirs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}_loh.csv"
    log:
        stdout = CFG["logs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purecn.stdout.log",
        stderr = CFG["logs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purecn.stderr.log"
    params:
        out_dir = CFG["dirs"]["purecn"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}",
        fun_segmentation = CFG["options"]["purecn_fun_segmentation"],
        min_purity = CFG["options"]["purecn_min_purity"],
        max_purity = CFG["options"]["purecn_max_purity"],
        max_copy_number = CFG["options"]["purecn_max_copy_number"],
        max_segments = CFG["options"]["purecn_max_segments"],
        # Genome arg: PureCN expects "hg19" / "hg38" labels; map grch37→hg19, grch38→hg38
        purecn_genome = lambda w: "hg19" if w.genome_build == "grch37" else ("hg38" if w.genome_build == "grch38" else w.genome_build)
    conda:
        CFG["conda_envs"]["purecn"]
    threads:
        CFG["threads"]["purecn"]
    resources:
        **CFG["resources"]["purecn"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.run_purecn}
            --seg_modelFinal {input.seg_modelFinal}
            --mutect_vcf     {input.mutect_vcf}
            --intervals      {input.intervals}
            --sample_id      {wildcards.tumour_id}
            --out_dir        {params.out_dir}
            --genome         {params.purecn_genome}
            --fun_segmentation {params.fun_segmentation}
            --min_purity     {params.min_purity}
            --max_purity     {params.max_purity}
            --max_copy_number {params.max_copy_number}
            --max_segments  {params.max_segments}
            > {log.stdout} 2> {log.stderr}
        """)


##### POSTPROCESSING #####


rule _gatk_cnv_seg_to_igv:
    input:
        called = str(rules._gatk_cnv_call_segments.output.called)
    output:
        igv = CFG["dirs"]["called"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.igv.seg"
    log:
        stderr = CFG["logs"]["called"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.seg_to_igv.stderr.log"
    params:
        sample = "{tumour_id}"
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["seg_to_igv"]
    resources:
        **CFG["resources"]["seg_to_igv"]
    shell:
        op.as_one_line("""
        (
          echo -e "ID\\tchrom\\tloc.start\\tloc.end\\tnum.mark\\tseg.mean";
          awk -v s="{params.sample}" 'BEGIN {{FS=OFS="\\t"}}
              !/^@/ && !/^CONTIG/ && NF >= 5 {{
                print s, $1, $2, $3, $4, $5;
              }}' {input.called}
        ) > {output.igv} 2> {log.stderr}
        """)


def _gatk_cnv_gene_matrix_inputs(wildcards):
    """Collect per-tumour called segs and PureCN gene CN files across all runs.
    Uses module-global _RUNS captured before op.cleanup_module pops it from CFG."""
    runs = _RUNS
    runs = runs[runs["tumour_genome_build"] == wildcards.genome_build]
    runs = runs[runs["tumour_seq_type"] == wildcards.seq_type]
    called = expand(
        str(rules._gatk_cnv_call_segments.output.called),
        zip,
        seq_type   = runs["tumour_seq_type"],
        genome_build = runs["tumour_genome_build"],
        tumour_id  = runs["tumour_sample_id"],
        normal_id  = runs["normal_sample_id"],
        pair_status = runs["pair_status"],
    )
    purecn_genes = expand(
        str(rules._gatk_cnv_purecn.output.genes_csv),
        zip,
        seq_type   = runs["tumour_seq_type"],
        genome_build = runs["tumour_genome_build"],
        tumour_id  = runs["tumour_sample_id"],
        normal_id  = runs["normal_sample_id"],
        pair_status = runs["pair_status"],
    )
    purecn_purity = expand(
        str(rules._gatk_cnv_purecn.output.purity_csv),
        zip,
        seq_type   = runs["tumour_seq_type"],
        genome_build = runs["tumour_genome_build"],
        tumour_id  = runs["tumour_sample_id"],
        normal_id  = runs["normal_sample_id"],
        pair_status = runs["pair_status"],
    )
    return {"called": called, "purecn_genes": purecn_genes, "purecn_purity": purecn_purity}


rule _gatk_cnv_gene_matrix:
    input:
        unpack(_gatk_cnv_gene_matrix_inputs),
        builder = ancient(CFG["inputs"]["build_gene_matrix"]),
        gtf     = reference_files("genomes/{genome_build}/annotations/gencode_annotation-" + CFG["options"]["gtf_release"] + ".gtf")
    output:
        relative = CFG["dirs"]["gene_matrix"] + "{seq_type}--{genome_build}/gene_relative_cn.tsv",
        absolute = CFG["dirs"]["gene_matrix"] + "{seq_type}--{genome_build}/gene_absolute_cn.tsv"
    log:
        stdout = CFG["logs"]["gene_matrix"] + "{seq_type}--{genome_build}/gene_matrix.stdout.log",
        stderr = CFG["logs"]["gene_matrix"] + "{seq_type}--{genome_build}/gene_matrix.stderr.log"
    conda:
        CFG["conda_envs"]["python"]
    threads:
        CFG["threads"]["gene_matrix"]
    resources:
        **CFG["resources"]["gene_matrix"]
    shell:
        op.as_one_line("""
        python {input.builder}
            --gtf {input.gtf}
            --called-segs {input.called}
            --purecn-genes {input.purecn_genes}
            --purity-csvs {input.purecn_purity}
            --relative-out {output.relative}
            --absolute-out {output.absolute}
            > {log.stdout} 2> {log.stderr}
        """)


def _gatk_cnv_summary_inputs(wildcards):
    """Uses module-global _RUNS captured before op.cleanup_module pops it from CFG."""
    runs = _RUNS
    runs = runs[runs["tumour_genome_build"] == wildcards.genome_build]
    runs = runs[runs["tumour_seq_type"] == wildcards.seq_type]
    return expand(
        str(rules._gatk_cnv_purecn.output.purity_csv),
        zip,
        seq_type   = runs["tumour_seq_type"],
        genome_build = runs["tumour_genome_build"],
        tumour_id  = runs["tumour_sample_id"],
        normal_id  = runs["normal_sample_id"],
        pair_status = runs["pair_status"],
    )


rule _gatk_cnv_purity_ploidy_summary:
    input:
        purity_csvs = _gatk_cnv_summary_inputs
    output:
        summary = CFG["dirs"]["gene_matrix"] + "{seq_type}--{genome_build}/purity_ploidy.tsv"
    log:
        stderr = CFG["logs"]["gene_matrix"] + "{seq_type}--{genome_build}/summary.stderr.log"
    threads: 1
    run:
        rows = []
        for csv_path in input.purity_csvs:
            try:
                df = pd.read_csv(csv_path)
                df["sample_id"] = os.path.splitext(os.path.basename(csv_path))[0]
                rows.append(df)
            except Exception as e:
                with open(log.stderr, "a") as fh:
                    fh.write(f"WARN: could not read {csv_path}: {e}\n")
        if rows:
            pd.concat(rows, ignore_index=True).to_csv(output.summary, sep="\t", index=False)
        else:
            # Emit empty file with header so downstream tools don't fail
            pd.DataFrame(columns=["sample_id", "Purity", "Ploidy"]).to_csv(output.summary, sep="\t", index=False)


##### OUTPUT SYMLINKS #####


rule _gatk_cnv_output_called_seg:
    input:
        called = str(rules._gatk_cnv_call_segments.output.called),
        igv    = str(rules._gatk_cnv_seg_to_igv.output.igv)
    output:
        called = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.called.seg",
        igv    = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.igv.seg"
    run:
        op.relative_symlink(input.called, output.called, in_module=True)
        op.relative_symlink(input.igv,    output.igv,    in_module=True)


rule _gatk_cnv_output_denoised:
    input:
        denoised_cr = str(rules._gatk_cnv_denoise_read_counts.output.denoised_cr)
    output:
        denoised_cr = CFG["dirs"]["outputs"] + "denoised/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoisedCR.tsv"
    run:
        op.relative_symlink(input.denoised_cr, output.denoised_cr, in_module=True)


rule _gatk_cnv_output_plots:
    input:
        denoised_plot = str(rules._gatk_cnv_plot_denoised.output.denoised_plot),
        modeled_plot  = str(rules._gatk_cnv_plot_modeled.output.modeled_plot)
    output:
        denoised_plot = CFG["dirs"]["outputs"] + "plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoised.png",
        modeled_plot  = CFG["dirs"]["outputs"] + "plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.modeled.png"
    run:
        op.relative_symlink(input.denoised_plot, output.denoised_plot, in_module=True)
        op.relative_symlink(input.modeled_plot,  output.modeled_plot,  in_module=True)


rule _gatk_cnv_output_purecn:
    input:
        purity_csv = str(rules._gatk_cnv_purecn.output.purity_csv),
        genes_csv  = str(rules._gatk_cnv_purecn.output.genes_csv),
        loh_csv    = str(rules._gatk_cnv_purecn.output.loh_csv)
    output:
        purity_csv = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.purity.csv",
        genes_csv  = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.genes.csv",
        loh_csv    = CFG["dirs"]["outputs"] + "purecn/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.loh.csv"
    run:
        op.relative_symlink(input.purity_csv, output.purity_csv, in_module=True)
        op.relative_symlink(input.genes_csv,  output.genes_csv,  in_module=True)
        op.relative_symlink(input.loh_csv,    output.loh_csv,    in_module=True)


rule _gatk_cnv_output_gene_matrix:
    input:
        relative = str(rules._gatk_cnv_gene_matrix.output.relative),
        absolute = str(rules._gatk_cnv_gene_matrix.output.absolute)
    output:
        relative = CFG["dirs"]["outputs"] + "gene_matrix/{seq_type}--{genome_build}/gene_relative_cn.tsv",
        absolute = CFG["dirs"]["outputs"] + "gene_matrix/{seq_type}--{genome_build}/gene_absolute_cn.tsv"
    run:
        op.relative_symlink(input.relative, output.relative, in_module=True)
        op.relative_symlink(input.absolute, output.absolute, in_module=True)


rule _gatk_cnv_output_summary:
    input:
        summary = str(rules._gatk_cnv_purity_ploidy_summary.output.summary)
    output:
        summary = CFG["dirs"]["outputs"] + "summary/{seq_type}--{genome_build}/purity_ploidy.tsv"
    run:
        op.relative_symlink(input.summary, output.summary, in_module=True)


##### AGGREGATOR #####


rule _gatk_cnv_all:
    input:
        # Per-tumour outputs
        expand(
            [
                str(rules._gatk_cnv_output_called_seg.output.called),
                str(rules._gatk_cnv_output_called_seg.output.igv),
                str(rules._gatk_cnv_output_denoised.output.denoised_cr),
                str(rules._gatk_cnv_output_plots.output.denoised_plot),
                str(rules._gatk_cnv_output_plots.output.modeled_plot),
                str(rules._gatk_cnv_output_purecn.output.purity_csv),
                str(rules._gatk_cnv_output_purecn.output.genes_csv),
                str(rules._gatk_cnv_output_purecn.output.loh_csv),
            ],
            zip,
            seq_type    = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id   = CFG["runs"]["tumour_sample_id"],
            normal_id   = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"],
        ),
        # Per-cohort (seq_type × genome_build) summary outputs
        expand(
            [
                str(rules._gatk_cnv_output_gene_matrix.output.relative),
                str(rules._gatk_cnv_output_gene_matrix.output.absolute),
                str(rules._gatk_cnv_output_summary.output.summary),
            ],
            zip,
            seq_type     = list(set(CFG["runs"]["tumour_seq_type"])),
            genome_build = list(set(CFG["runs"]["tumour_genome_build"])),
        )


##### CLEANUP #####


op.cleanup_module(CFG)
