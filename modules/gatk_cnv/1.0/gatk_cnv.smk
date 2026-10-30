#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import os
import sys
import pandas as pd

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
# `CFG` is a shortcut to `config["lcr-modules"]["gatk_cnv"]`
CFG = op.setup_module(
    name = "gatk_cnv",
    version = "1.0",
    subdirectories = ["inputs", "intervals", "read_counts", "pon", "allelic_counts", "denoise", "model_segments",
                      "seg", "convert_coordinates", "fill_regions", "normalize", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gatk_cnv_input_bam,
    _gatk_cnv_input_references,
    _gatk_cnv_filter_dropouts,
    _gatk_cnv_to_seg,
    _gatk_cnv_normalize_projection,
    _gatk_cnv_output_projection,
    _gatk_cnv_output,
    _gatk_cnv_output_plot,
    _gatk_cnv_all


# Java heap for GATK: most of the memory requested for the job
def _gatk_cnv_java_mem(step):
    return str(int(CFG["resources"][step]["mem_mb"] * 0.8)) + "m"


# ModelSegments smoothing/segmentation arguments shared by both pairing modes
_GATK_CNV_MS_OPTS = op.as_one_line("""
    --number-of-changepoints-penalty-factor {changepoints_penalty_factor}
    --number-of-smoothing-iterations-per-fit {smoothing_iterations_per_fit}
    --smoothing-credible-interval-threshold-copy-ratio {smoothing_credible_interval_threshold}
    --smoothing-credible-interval-threshold-allele-fraction {smoothing_credible_interval_threshold}
    {opts}
""".format(**CFG["options"]["model_segments"]))


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gatk_cnv_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


# Reference genome (CRAM decoding and every GATK step use the same FASTA)
rule _gatk_cnv_input_references:
    input:
        fasta = ancient(reference_files("genomes/{genome_build}/genome_fasta/genome.fa")),
        fai = ancient(reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")),
        gatk_dict = ancient(reference_files("genomes/{genome_build}/genome_fasta/genome.dict"))
    output:
        fasta = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa",
        fai = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.fa.fai",
        gatk_dict = CFG["dirs"]["inputs"] + "references/{genome_build}/genome.dict"
    run:
        op.absolute_symlink(input.fasta, output.fasta)
        op.absolute_symlink(input.fai, output.fai)
        op.absolute_symlink(input.gatk_dict, output.gatk_dict)


# Padded capture targets (bin length 0: one interval per target), restricted to
# the contigs matching `chromosome_regex`
rule _gatk_cnv_preprocess_intervals:
    input:
        baits = ancient(reference_files("genomes/{genome_build}/capture_space/{capture_space}.baits.bed")),
        fasta = str(rules._gatk_cnv_input_references.output.fasta),
        gatk_dict = str(rules._gatk_cnv_input_references.output.gatk_dict)
    output:
        all_contigs = temp(CFG["dirs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/targets.all_contigs.interval_list"),
        intervals = CFG["dirs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/targets.preprocessed.interval_list"
    log:
        stderr = CFG["logs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/preprocess_intervals.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("preprocess"),
        padding = CFG["options"]["preprocess"]["padding"],
        chromosome_regex = CFG["options"]["preprocess"]["chromosome_regex"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["preprocess"]
    resources:
        **CFG["resources"]["preprocess"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" PreprocessIntervals
            -R {input.fasta} -L {input.baits} --bin-length 0 --padding {params.padding}
            --interval-merging-rule OVERLAPPING_ONLY -O {output.all_contigs} 2> {log.stderr}
        &&
        awk -v re='{params.chromosome_regex}' '/^@/ {{print; next}} $1 ~ re {{print}}'
            {output.all_contigs} > {output.intervals}
        """)


# GC content per target, for the panel of normals' GC correction
rule _gatk_cnv_annotate_intervals:
    input:
        intervals = str(rules._gatk_cnv_preprocess_intervals.output.intervals),
        fasta = str(rules._gatk_cnv_input_references.output.fasta)
    output:
        annotated = CFG["dirs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/targets.annotated.tsv"
    log:
        stderr = CFG["logs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/annotate_intervals.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("preprocess")
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["preprocess"]
    resources:
        **CFG["resources"]["preprocess"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" AnnotateIntervals
            -R {input.fasta} -L {input.intervals} --interval-merging-rule OVERLAPPING_ONLY
            -O {output.annotated} 2> {log.stderr}
        """)


# Common biallelic SNP sites (gnomAD AF > min_snp_af) within the targets, where
# allele counts are collected. -R supplies the sequence dictionary, since the
# gnomAD VCF header may lack ##contig lines.
rule _gatk_cnv_common_snps:
    input:
        intervals = str(rules._gatk_cnv_preprocess_intervals.output.intervals),
        fasta = str(rules._gatk_cnv_input_references.output.fasta),
        gnomad = ancient(reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz"))
    output:
        vcf = CFG["dirs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/common_snps.vcf.gz",
        tbi = CFG["dirs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/common_snps.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["intervals"] + "{seq_type}--{genome_build}/{capture_space}/common_snps.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("preprocess"),
        min_af = CFG["options"]["min_snp_af"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["preprocess"]
    resources:
        **CFG["resources"]["preprocess"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" SelectVariants
            -R {input.fasta} -V {input.gnomad} -L {input.intervals}
            --select-type-to-include SNP --restrict-alleles-to BIALLELIC
            -select "AF > {params.min_af}" -O {output.vcf} 2> {log.stderr}
        """)


# Read counts per target, for tumours and panel-of-normals normals alike
rule _gatk_cnv_collect_read_counts:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        intervals = str(rules._gatk_cnv_preprocess_intervals.output.intervals),
        fasta = str(rules._gatk_cnv_input_references.output.fasta)
    output:
        counts = CFG["dirs"]["read_counts"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.counts.hdf5"
    log:
        stderr = CFG["logs"]["read_counts"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("read_counts")
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["read_counts"]
    resources:
        **CFG["resources"]["read_counts"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" CollectReadCounts
            -I {input.bam} -R {input.fasta} -L {input.intervals}
            --interval-merging-rule OVERLAPPING_ONLY -O {output.counts} 2> {log.stderr}
        """)


# Normals in the samples table for this seq_type, genome_build and capture_space,
# minus `pon_exclude`. Tumours' own normals are included (one panel per capture space).
def _gatk_cnv_pon_normals(wildcards):
    CFG = config["lcr-modules"]["gatk_cnv"]
    tbl = CFG["samples"]
    normals = tbl[
        (tbl.tissue_status.str.lower() == "normal") &
        (tbl.seq_type == wildcards.seq_type) &
        (tbl.genome_build == wildcards.genome_build) &
        (tbl.capture_space == wildcards.capture_space) &
        (~tbl.sample_id.isin(CFG["options"]["pon_exclude"]))
    ]["sample_id"].drop_duplicates().tolist()
    if len(normals) < CFG["options"]["pon_min_normals"]:
        raise ValueError(
            f"Only {len(normals)} normals for {wildcards.seq_type}--{wildcards.genome_build}--{wildcards.capture_space}; "
            f"the panel of normals needs at least {CFG['options']['pon_min_normals']} (options: pon_min_normals)."
        )
    return expand(
        str(rules._gatk_cnv_collect_read_counts.output.counts),
        sample_id = sorted(normals),
        allow_missing = True
    )


rule _gatk_cnv_create_pon:
    input:
        counts = _gatk_cnv_pon_normals,
        annotated = str(rules._gatk_cnv_annotate_intervals.output.annotated)
    output:
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/pon.hdf5",
        normals = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/pon_normals.txt"
    log:
        stderr = CFG["logs"]["pon"] + "{seq_type}--{genome_build}/{capture_space}/create_pon.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("pon"),
        inputs = lambda w, input: " ".join("-I " + f for f in input.counts)
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["pon"]
    resources:
        **CFG["resources"]["pon"]
    shell:
        op.as_one_line("""
        for f in {input.counts}; do basename $f .counts.hdf5; done > {output.normals}
        &&
        gatk --java-options "-Xmx{params.java_mem}" CreateReadCountPanelOfNormals
            {params.inputs} --annotated-intervals {input.annotated}
            -O {output.pon} 2> {log.stderr}
        """)


# Ref/alt counts at the common SNP sites (tumours and matched normals)
rule _gatk_cnv_collect_allelic_counts:
    input:
        bam = str(rules._gatk_cnv_input_bam.output.bam),
        bai = str(rules._gatk_cnv_input_bam.output.bai),
        snps = str(rules._gatk_cnv_common_snps.output.vcf),
        fasta = str(rules._gatk_cnv_input_references.output.fasta)
    output:
        counts = CFG["dirs"]["allelic_counts"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.allelicCounts.tsv"
    log:
        stderr = CFG["logs"]["allelic_counts"] + "{seq_type}--{genome_build}/{capture_space}/{sample_id}.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("allelic_counts")
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["allelic_counts"]
    resources:
        **CFG["resources"]["allelic_counts"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" CollectAllelicCounts
            -I {input.bam} -R {input.fasta} -L {input.snps} -O {output.counts} 2> {log.stderr}
        """)


# Standardize and denoise tumour read counts against the panel of normals
rule _gatk_cnv_denoise:
    input:
        counts = CFG["dirs"]["read_counts"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.counts.hdf5",
        pon = str(rules._gatk_cnv_create_pon.output.pon)
    output:
        standardized = CFG["dirs"]["denoise"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.standardizedCR.tsv",
        denoised = CFG["dirs"]["denoise"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.denoisedCR.tsv"
    log:
        stderr = CFG["logs"]["denoise"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("denoise"),
        opts = CFG["options"]["denoise"]["opts"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["denoise"]
    resources:
        **CFG["resources"]["denoise"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" DenoiseReadCounts
            -I {input.counts} --count-panel-of-normals {input.pon} {params.opts}
            --standardized-copy-ratios {output.standardized}
            --denoised-copy-ratios {output.denoised} 2> {log.stderr}
        """)


# Drop outlier targets before segmentation: zero-coverage dropouts (denoised
# log2 below min_log2_copy_ratio) and high-side artefacts (above
# max_log2_copy_ratio), which otherwise each become a few-target segment
rule _gatk_cnv_filter_dropouts:
    input:
        denoised = str(rules._gatk_cnv_denoise.output.denoised)
    output:
        filtered = CFG["dirs"]["denoise"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.denoisedCR.filtered.tsv"
    params:
        min_log2 = CFG["options"]["min_log2_copy_ratio"],
        max_log2 = CFG["options"]["max_log2_copy_ratio"]
    shell:
        op.as_one_line("""
        awk -F'\\t' -v min={params.min_log2} -v max={params.max_log2}
            '/^@/ || $1 == "CONTIG" || ($4 >= min && $4 <= max)' {input.denoised} > {output.filtered}
        """)


# Joint copy-ratio and allele-fraction segmentation; heterozygous sites are
# genotyped from the matched normal
rule _gatk_cnv_model_segments:
    input:
        denoised = str(rules._gatk_cnv_filter_dropouts.output.filtered),
        tumour_counts = CFG["dirs"]["allelic_counts"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.allelicCounts.tsv",
        normal_counts = CFG["dirs"]["allelic_counts"] + "{seq_type}--{genome_build}/{capture_space}/{normal_id}.allelicCounts.tsv"
    output:
        model_final = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.seg",
        cr_seg = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cr.seg",
        hets = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.hets.tsv"
    log:
        stderr = CFG["logs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/model_segments.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("model_segments"),
        outdir = lambda w, output: os.path.dirname(output.model_final),
        min_case = CFG["options"]["model_segments"]["min_total_allele_count_case_paired"],
        opts = _GATK_CNV_MS_OPTS
    wildcard_constraints:
        pair_status = "matched|unmatched"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["model_segments"]
    resources:
        **CFG["resources"]["model_segments"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" ModelSegments
            --denoised-copy-ratios {input.denoised}
            --allelic-counts {input.tumour_counts} --normal-allelic-counts {input.normal_counts}
            --minimum-total-allele-count-case {params.min_case}
            {params.opts}
            --output {params.outdir} --output-prefix {wildcards.tumour_id} 2> {log.stderr}
        """)


# Tumour-only: heterozygous sites are genotyped from the tumour itself, which
# needs read depth. Copy ratios are as reliable as with a normal; allele
# fractions are not, in LOH regions of high-purity tumours.
rule _gatk_cnv_model_segments_no_normal:
    input:
        denoised = str(rules._gatk_cnv_filter_dropouts.output.filtered),
        tumour_counts = CFG["dirs"]["allelic_counts"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.allelicCounts.tsv"
    output:
        model_final = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.seg",
        cr_seg = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cr.seg",
        hets = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.hets.tsv"
    log:
        stderr = CFG["logs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/model_segments.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("model_segments"),
        outdir = lambda w, output: os.path.dirname(output.model_final),
        min_case = CFG["options"]["model_segments"]["min_total_allele_count_case_no_normal"],
        opts = _GATK_CNV_MS_OPTS
    wildcard_constraints:
        pair_status = "no_normal"
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["model_segments"]
    resources:
        **CFG["resources"]["model_segments"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" ModelSegments
            --denoised-copy-ratios {input.denoised}
            --allelic-counts {input.tumour_counts}
            --minimum-total-allele-count-case {params.min_case}
            {params.opts}
            --output {params.outdir} --output-prefix {wildcards.tumour_id} 2> {log.stderr}
        """)


# Neutral / amplified / deleted calls for the copy-ratio segments
rule _gatk_cnv_call_segments:
    input:
        cr_seg = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.cr.seg"
    output:
        called = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.called.seg"
    log:
        stderr = CFG["logs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/call_segments.stderr.log"
    params:
        java_mem = _gatk_cnv_java_mem("post_process"),
        opts = CFG["options"]["call"]["opts"]
    conda:
        CFG["conda_envs"]["gatk"]
    threads:
        CFG["threads"]["post_process"]
    resources:
        **CFG["resources"]["post_process"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.java_mem}" CallCopyRatioSegments
            -I {input.cr_seg} {params.opts} -O {output.called} 2> {log.stderr}
        """)


# Genome-wide plot: denoised copy ratio and allele fraction with the final
# segments; targets removed before segmentation are shown in grey
rule _gatk_cnv_plot:
    input:
        denoised = str(rules._gatk_cnv_denoise.output.denoised),
        hets = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.hets.tsv",
        model_final = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.seg"
    output:
        png = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.png"
    log:
        stderr = CFG["logs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/plot.stderr.log"
    params:
        script = CFG["software"]["plot_script"],
        min_log2 = CFG["options"]["min_log2_copy_ratio"],
        max_log2 = CFG["options"]["max_log2_copy_ratio"],
        chromosome_regex = CFG["options"]["preprocess"]["chromosome_regex"]
    conda:
        CFG["conda_envs"]["plot"]
    threads:
        CFG["threads"]["post_process"]
    resources:
        **CFG["resources"]["post_process"]
    shell:
        op.as_one_line("""
        Rscript {params.script} --denoised {input.denoised} --hets {input.hets}
            --segments {input.model_final} --sample {wildcards.tumour_id}
            --min_log2 {params.min_log2} --max_log2 {params.max_log2}
            --chromosome_regex '{params.chromosome_regex}' --out {output.png} 2> {log.stderr}
        """)


# Standardized seg (ID chrom start end module LOH_flag CN log.ratio), the
# layout fill_segments, cnv_master and GAMBLR expect. GATK segments are
# 1-based closed; starts are converted to 0-based. CN = round(2 * 2^log2)
# (not purity-corrected); LOH_flag is NA.
rule _gatk_cnv_to_seg:
    input:
        model_final = CFG["dirs"]["model_segments"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}/{tumour_id}.modelFinal.seg"
    output:
        seg = CFG["dirs"]["seg"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}.seg"
    params:
        tool = CFG["options"]["tool_name"]
    run:
        segs = pd.read_csv(input.model_final, sep = "\t", comment = "@")
        segs = segs[segs["LOG2_COPY_RATIO_POSTERIOR_50"].notna()]
        log2 = segs["LOG2_COPY_RATIO_POSTERIOR_50"]
        out = pd.DataFrame({
            "ID": wildcards.tumour_id,
            "chrom": segs["CONTIG"].astype(str),
            "start": segs["START"] - 1,
            "end": segs["END"],
            "module": params.tool,
            "LOH_flag": "NA",
            "CN": (2 * 2 ** log2).round().astype(int),
            "log.ratio": log2.round(4)
        })
        out.to_csv(output.seg, sep = "\t", index = False, na_rep = "NA")


def _gatk_cnv_get_chain(wildcards):
    if "38" in str(wildcards.genome_build):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")


# Convert the coordinates of the seg file to the other genome build
rule _gatk_cnv_convert_coordinates:
    input:
        native = str(rules._gatk_cnv_to_seg.output.seg),
        chain = _gatk_cnv_get_chain
    output:
        lifted = CFG["dirs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}--{normal_id}--{pair_status}.lifted_{chain}.stderr.log"
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    threads:
        CFG["threads"]["post_process"]
    resources:
        **CFG["resources"]["post_process"]
    shell:
        op.as_one_line("""
        bash {params.liftover_script}
        SEG
        {input.native}
        {output.lifted}
        {input.chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)


# Look up a tumour's genome build and capture space from the runs table
def _gatk_cnv_run_attr(wildcards, attr):
    tbl = config["lcr-modules"]["gatk_cnv"]["runs"]
    return tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_" + attr].tolist()[0]


def _gatk_cnv_prepare_projection(wildcards):
    this_build = _gatk_cnv_run_attr(wildcards, "genome_build")
    this_space = _gatk_cnv_run_attr(wildcards, "capture_space")
    native = str(rules._gatk_cnv_to_seg.output.seg).replace("{genome_build}", this_build).replace("{capture_space}", this_space)
    lifted = str(rules._gatk_cnv_convert_coordinates.output.lifted).replace("{genome_build}", this_build).replace("{capture_space}", this_space)
    if "38" in this_build:
        return {"hg38_projection": native, "grch37_projection": lifted.replace("{chain}", "hg38ToHg19")}
    else:
        return {"grch37_projection": native, "hg38_projection": lifted.replace("{chain}", "hg19ToHg38")}


# Fill the missing segments of seg files with neutral regions to complete the genome coverage
rule _gatk_cnv_fill_segments:
    input:
        unpack(_gatk_cnv_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}_fill_segments.stderr.log"
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    threads:
        CFG["threads"]["post_process"]
    resources:
        **CFG["resources"]["post_process"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id}--{wildcards.normal_id} on $(hostname) at $(date)" > {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.hg38.bed
        {input.hg38_projection}
        {params.path}src/blacklisted.hg38.bed
        {output.hg38_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        """)


# Normalize chr prefix of the output file
rule _gatk_cnv_normalize_projection:
    input:
        filled = CFG["dirs"]["fill_regions"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg",
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    wildcard_constraints:
        projection = "|".join(CFG["output"]["requested_projections"])
    run:
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names = ["chromosome"], header = None)
        seg = pd.read_csv(input.filled, sep = "\t")
        chrom = seg["chrom"].astype(str).str.replace("chr", "")
        if "chr" in str(chromosomes["chromosome"][0]):
            chrom = "chr" + chrom
        seg["chrom"] = chrom
        seg.to_csv(output.projection, sep = "\t", index = False, na_rep = "NA")


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gatk_cnv_output_projection:
    input:
        projection = str(rules._gatk_cnv_normalize_projection.output.projection)
    output:
        projection = CFG["dirs"]["outputs"] + "seg/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{tool}.{projection}.seg"
    wildcard_constraints:
        projection = "|".join(CFG["output"]["requested_projections"])
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)


# Outputs without the capture_space wildcard, for consistency with other CNV modules
def _gatk_cnv_drop_capture_space_wc(wildcards):
    this_space = _gatk_cnv_run_attr(wildcards, "capture_space")
    files = {
        "seg": str(rules._gatk_cnv_to_seg.output.seg),
        "model_final": str(rules._gatk_cnv_model_segments.output.model_final),
        "called": str(rules._gatk_cnv_call_segments.output.called),
        "hets": str(rules._gatk_cnv_model_segments.output.hets),
        "denoised": str(rules._gatk_cnv_denoise.output.denoised)
    }
    return {k: v.replace("{capture_space}", this_space) for k, v in files.items()}


rule _gatk_cnv_output:
    input:
        unpack(_gatk_cnv_drop_capture_space_wc)
    output:
        seg = CFG["dirs"]["outputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.seg",
        model_final = CFG["dirs"]["outputs"] + "model_segments/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.modelFinal.seg",
        called = CFG["dirs"]["outputs"] + "model_segments/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.called.seg",
        hets = CFG["dirs"]["outputs"] + "model_segments/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.hets.tsv",
        denoised = CFG["dirs"]["outputs"] + "denoised/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.denoisedCR.tsv"
    run:
        for key in ["seg", "model_final", "called", "hets", "denoised"]:
            op.relative_symlink(input[key], output[key], in_module = True)


def _gatk_cnv_plot_wc(wildcards):
    this_space = _gatk_cnv_run_attr(wildcards, "capture_space")
    return {"png": str(rules._gatk_cnv_plot.output.png).replace("{capture_space}", this_space)}


rule _gatk_cnv_output_plot:
    input:
        unpack(_gatk_cnv_plot_wc)
    output:
        png = CFG["dirs"]["outputs"] + "plots/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.modelFinal.png"
    run:
        op.relative_symlink(input.png, output.png, in_module = True)


# Plots are optional (output: plot)
_GATK_CNV_PLOTS = expand(
    str(rules._gatk_cnv_output_plot.output.png),
    zip,
    seq_type = CFG["runs"]["tumour_seq_type"],
    genome_build = CFG["runs"]["tumour_genome_build"],
    tumour_id = CFG["runs"]["tumour_sample_id"],
    normal_id = CFG["runs"]["normal_sample_id"],
    pair_status = CFG["runs"]["pair_status"]
) if CFG["output"]["plot"] else []


# Generates the target sentinels for each run, which generate the symlinks
rule _gatk_cnv_all:
    input:
        _GATK_CNV_PLOTS,
        expand(
            [
                str(rules._gatk_cnv_output.output.seg),
                str(rules._gatk_cnv_output.output.model_final),
                str(rules._gatk_cnv_output.output.called),
                str(rules._gatk_cnv_output.output.hets),
                str(rules._gatk_cnv_output.output.denoised)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type = CFG["runs"]["tumour_seq_type"],
            genome_build = CFG["runs"]["tumour_genome_build"],
            tumour_id = CFG["runs"]["tumour_sample_id"],
            normal_id = CFG["runs"]["normal_sample_id"],
            pair_status = CFG["runs"]["pair_status"]
        ),
        expand(
            expand(
                str(rules._gatk_cnv_output_projection.output.projection),
                zip,  # Run expand() with zip(), not product()
                seq_type = CFG["runs"]["tumour_seq_type"],
                tumour_id = CFG["runs"]["tumour_sample_id"],
                normal_id = CFG["runs"]["normal_sample_id"],
                pair_status = CFG["runs"]["pair_status"],
                allow_missing = True
            ),
            tool = CFG["options"]["tool_name"],
            projection = CFG["output"]["requested_projections"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
