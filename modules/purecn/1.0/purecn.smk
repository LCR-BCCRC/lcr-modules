#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Jasper Wong
# Module Author:    Jasper Wong
# Contributors:     Sierra Gillis


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import glob
import os
import re
import sys
from pathlib import Path

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
# `CFG` is a shortcut to `config["lcr-modules"]["purecn"]`
CFG = op.setup_module(
    name = "purecn",
    version = "1.0",
    subdirectories = ["inputs", "mutect2", "coverage", "pureCN_cnvkit", "pureCN_denovo", "pureCN_final_result", "cnv2igv", "convert_coordinates", "fill_regions", "normalize", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _purecn_input_bam,
    _purecn_symlink_cnvkit_seg,
    _purecn_symlink_cnvkit_cnr,
    _purecn_setup_blacklist_bed,
    _purecn_input_chroms_withY,
    _purecn_symlink_fasta,
    _purecn_symlink_intervals,
    _purecn_symlink_gatk_intervals,
    _purecn_gatk_interval_list_chrom,
    _purecn_symlink_gatk_intervals_targets,
    _purecn_symlink_mutect2_pon,
    _purecn_symlink_mutect2_pon_tbi,
    _purecn_symlink_database_cnvkit,
    _purecn_symlink_database_denovo,
    _purecn_make_outdirs_cnvkit_mode,
    _purecn_symlink_cnvkit_mode_run_result,
    _purecn_make_outdirs_denovo_mode,
    _purecn_symlink_denovo_mode_run_result,
    _purecn_record_seg_fn,
    _purecn_output_projection,
    _purecn_output,
    _purecn_best_seg,
    _purecn_best_seg_info,
    _purecn_all


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

rule _purecn_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)

# Import CNVkit seg and coverage
rule _purecn_symlink_cnvkit_seg:
    input:
        seg = CFG["inputs"]["cnvkit_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
    run:
        op.absolute_symlink(input.seg, output.seg)

rule _purecn_symlink_cnvkit_cnr:
    input:
        cnr = CFG["inputs"]["cnvkit_cnr"]
    output:
        cnr = CFG["dirs"]["inputs"] + "cnr/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.cnr"
    run:
        op.absolute_symlink(input.cnr, output.cnr)

rule _purecn_setup_blacklist_bed:
    input:
        blacklist = ancient(reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed"))
    output:
        blacklist = CFG["dirs"]["inputs"] + "references/{genome_build}/repeatmasker.{genome_build}.bed"
    shell:
        """
            awk '{{print $1"\t"$2"\t"$3}}' {input.blacklist} > {output.blacklist}
        """

checkpoint _purecn_input_chroms_withY:
    input:
        txt = ancient(reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes_withY.txt"))
    output:
        txt = CFG["dirs"]["inputs"] + "references/{genome_build}/main_chromosomes_withY.txt"
    run:
        op.absolute_symlink(input.txt, output.txt)

# Symlink fasta and it's indexes (for compatibility with PURECN Rscripts and MuTect2/GATK)
rule _purecn_symlink_fasta:
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

# Intervals made by PureCN in _panel_of_normals
rule _purecn_symlink_intervals:
    input:
        intervals = ancient(CFG["inputs"]["purecn_intervals"])
    output:
        intervals = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals.txt"
    run:
        op.absolute_symlink(input.intervals, output.intervals)

# Used for calculating coverage with GATK since PureCN cannot take crams directly
rule _purecn_symlink_gatk_intervals:
    input:
        gatk_intervals = ancient(CFG["inputs"]["gatk_intervals"])
    output:
         gatk_intervals = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk.list"
    run:
        op.absolute_symlink(input.gatk_intervals, output.gatk_intervals)

# Split by chrom, used later in GATK depthOfCoverage
rule _purecn_gatk_interval_list_chrom:
    input:
        gatk_intervals = str(rules._purecn_symlink_gatk_intervals.output.gatk_intervals)
    output:
        chrom_int = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_{chrom}.intervals_gatk.list"
    log:
        CFG["logs"]["inputs"] + "{seq_type}--{genome_build}/{capture_space}/purecn_gatk_intervals_{chrom}.log"
    shell:
        op.as_one_line(
        """
            num_intervals=$( {{ egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} || true; }} | wc -l );
            if [[ $num_intervals -eq 0 ]]; then
                echo "No intervals found for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log};
                echo "{wildcards.chrom}" > {output.chrom_int};
            else
                echo "Found $num_intervals intervals for chromosome {wildcards.chrom} in {input.gatk_intervals}" | tee {log};
                egrep -i '^{wildcards.chrom}:.*-.*' {input.gatk_intervals} > {output.chrom_int};
            fi
        """
        )

# The sorted one with off-target regions dropped for MuTect2
rule _purecn_symlink_gatk_intervals_targets:
    input:
        gatk_targets = ancient(CFG["inputs"]["gatk_targets"])
    output:
        gatk_targets = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/baits_{capture_space}_intervals_gatk_targets_sorted.list"
    run:
        op.absolute_symlink(input.gatk_targets, output.gatk_targets)

# Symlink pon vcf made for running Mutect2
rule _purecn_symlink_mutect2_pon:
    input:
        mutect2_pon = ancient(CFG["inputs"]["mutect2_pon"])
    output:
        mutect2_pon = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/{capture_space}_mutect2_pon.vcf.gz"
    run:
        op.absolute_symlink(input.mutect2_pon, output.mutect2_pon)

# Symlink pon vcf index made for running Mutect2
rule _purecn_symlink_mutect2_pon_tbi:
    input:
        mutect2_pon_tbi = ancient(CFG["inputs"]["mutect2_pon_tbi"])
    output:
        mutect2_pon_tbi = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/{capture_space}_mutect2_pon.vcf.gz.tbi"
    run:
        op.absolute_symlink(input.mutect2_pon_tbi, output.mutect2_pon_tbi)

# Symlink mapping_bias file used in cnvkit mode
rule _purecn_symlink_database_cnvkit:
    input:
        mapping_bias = ancient(CFG["inputs"]["cnvkit_mapping_bias"])
    output:
        mapping_bias = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/purecn_cnvkit_normal/mapping_bias_{genome_build}_{capture_space}.rds"
    run:
        op.absolute_symlink(input.mapping_bias, output.mapping_bias)

# Symlink mapping bias and normalDB used in denovo mode
rule _purecn_symlink_database_denovo:
    input:
        mapping_bias = ancient(CFG["inputs"]["denovo_mapping_bias"]),
        normal_db = ancient(CFG["inputs"]["denovo_normal_db"])
    output:
        mapping_bias = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/mapping_bias_{genome_build}_{capture_space}.rds",
        normal_db = CFG["dirs"]["inputs"] + "panel_of_normals/{seq_type}--{genome_build}/{capture_space}/purecn_denovo_normal/normalDB_{genome_build}_{capture_space}.rds"
    run:
        op.absolute_symlink(input.mapping_bias, output.mapping_bias)
        op.absolute_symlink(input.normal_db, output.normal_db)

# Run Mutect2
rule _purecn_mutect2:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        dbsnp = ancient(reference_files("genomes/{genome_build}/variation/dbsnp.common_all-151.vcf.gz")),
        fasta = str(rules._purecn_symlink_fasta.output.fasta),
        gatk_dict = str(rules._purecn_symlink_fasta.output.gatk_dict),
        gnomad = ancient(reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")),
        pon = str(rules._purecn_symlink_mutect2_pon.output.mutect2_pon),
        pon_tbi = str(rules._purecn_symlink_mutect2_pon_tbi.output.mutect2_pon_tbi),
        target_regions = str(rules._purecn_symlink_gatk_intervals_targets.output.gatk_targets)
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz.tbi"),
        stats = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.vcf.gz.stats"),
        f1r2 = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.f1r2.tar.gz")
    resources:
        **CFG["resources"]["mutect"]
    threads: 1 # MuTect2 doesn't support multi-threaded and adding multiple CPUs can lead to memory bloat
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2"]["mutect2_opts"],
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/mutect2_{chrom}.log"
    conda:
        CFG["conda_envs"]["mutect"]
    shell:
        op.as_one_line("""
            if [[ $(egrep "^{wildcards.chrom}:" {input.target_regions} | wc -l) -eq 0 ]]; then
                echo "No intervals found for chromosome {wildcards.chrom} in {input.target_regions}" | tee {log};
                gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
                    Mutect2 {params.opts}
                    --genotype-germline-sites true
                    --genotype-pon-sites true
                    --germline-resource {input.gnomad}
                    -R {input.fasta}
                    -L {wildcards.chrom}:1-100
                    -pon {input.pon}
                    -I {input.bam}
                    -O {output.vcf}
                    --f1r2-tar-gz {output.f1r2}
                    >> {log} 2>&1;
            else
                echo "Found intervals for chromosome {wildcards.chrom} in {input.target_regions}" | tee {log};
                gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
                    Mutect2 {params.opts}
                    --genotype-germline-sites true
                    --genotype-pon-sites true
                    --germline-resource {input.gnomad}
                    -R {input.fasta}
                    -L {wildcards.chrom}
                    -L {input.target_regions}
                    -isr INTERSECTION
                    -pon {input.pon}
                    -I {input.bam}
                    -O {output.vcf}
                    --f1r2-tar-gz {output.f1r2}
                    >> {log} 2>&1;
            fi
        """)

# Concat vcfs and indexes across chroms per sample
def _get_mutect2_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    vcfs = expand(
        str(rules._purecn_mutect2.output.vcf),
        chrom = chrs,
        allow_missing = True
    )
    return(vcfs)

def _get_mutect2_chr_tbis(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    tbis = expand(
        str(rules._purecn_mutect2.output.tbi),
        chrom = chrs,
        allow_missing = True
    )
    return(tbis)

rule _purecn_mutect2_concat_vcf_per_sample:
    input:
        vcf = _get_mutect2_chr_vcfs,
        tbi = _get_mutect2_chr_tbis,
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_tmp.vcf.gz"),
        tbi = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_tmp.vcf.gz.tbi")
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_concat_vcf.log"
    resources:
        **CFG["resources"]["post_vcf"]
    threads: 1
    conda:
        CFG["conda_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools concat {input.vcf} -Oz -o {output.vcf} &> {log}  &&
         tabix -p vcf {output.vcf} &>> {log}
        """)

def _get_mutect2_chr_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    stats = expand(
        str(rules._purecn_mutect2.output.stats),
        chrom = chrs,
        allow_missing = True
    )
    return(stats)

# Merge mutect2 stats per chrom, for FilterMutectCalls rule
rule _purecn_mutect2_merge_stats_per_sample:
    input:
        stats = _get_mutect2_chr_stats
    output:
        stats = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_tmp.vcf.gz.stats"
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_merge_stats.log"
    params:
        mem_mb =  lambda wildcards, resources: int(resources.mem_mb * 0.8),
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["post_vcf"]
    threads: 1
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}" MergeMutectStats $(for i in {input.stats}; do echo -n "-stats $i "; done)
        -O {output.stats} > {log} 2>&1
        """)

# Get pileup summaries, for GATK CalculateContamination
rule _purecn_pileup_summaries:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        snps = ancient(reference_files("genomes/{genome_build}/gatk/mutect2_small_exac.{genome_build}.vcf.gz")),
        fasta = str(rules._purecn_symlink_fasta.output.fasta),
        gatk_dict = str(rules._purecn_symlink_fasta.output.gatk_dict)
    output:
        pileup = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/pileupSummary.table"
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_pileupSummary.log"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["post_vcf"]
    threads:
        CFG["threads"]["post_vcf"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
            GetPileupSummaries
            -I {input.bam}
            -R {input.fasta}
            -V {input.snps}
            -L {input.snps}
            -O {output.pileup}
            > {log} 2>&1
        """)

# Calculate contamination
rule _purecn_calc_contamination:
    input:
        pileup = str(rules._purecn_pileup_summaries.output.pileup)
    output:
        segments =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/segments.table",
        contamination =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/contamination.table"
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_calc_contam.log"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["post_vcf"]
    threads:
        CFG["threads"]["post_vcf"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
            CalculateContamination
            -I {input.pileup}
            -tumor-segmentation {output.segments}
            -O {output.contamination}
            > {log} 2>&1
        """)

# Learn read orientation model
def _get_mutect2_chr_f1r2(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    f1r2 = expand(
        str(rules._purecn_mutect2.output.f1r2),
        chrom = chrs,
        allow_missing = True
    )
    return(f1r2)

rule _purecn_learn_orient_model:
    input:
        f1r2 = _get_mutect2_chr_f1r2
    output:
        model =  CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/read-orientation-model.tar.gz"
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_learn_orient_model.log"
    conda:
        CFG["conda_envs"]["mutect"]
    resources:
        **CFG["resources"]["post_vcf"]
    threads:
        CFG["threads"]["post_vcf"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    shell:
        op.as_one_line("""
        inputs=$(for input in {input.f1r2}; do printf -- "-I $input "; done);
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
        LearnReadOrientationModel
        $inputs -O {output.model}
        > {log} 2>&1
        """)

# Marks variants filtered or PASS annotations
rule _purecn_annotate_vcf:
    input:
        vcf = str(rules._purecn_mutect2_concat_vcf_per_sample.output.vcf),
        tbi = str(rules._purecn_mutect2_concat_vcf_per_sample.output.tbi),
        stat = str(rules._purecn_mutect2_merge_stats_per_sample.output.stats),
        segments = str(rules._purecn_calc_contamination.output.segments),
        contamination = str(rules._purecn_calc_contamination.output.contamination),
        model = str(rules._purecn_learn_orient_model.output.model),
        fasta = str(rules._purecn_symlink_fasta.output.fasta)
    output:
        vcf = temp(CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/annotated.vcf.gz")
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_annotate_vcf.log",
    resources:
        **CFG["resources"]["post_vcf"]
    threads: 1
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["mutect2"]["annotate"],
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    conda:
        CFG["conda_envs"]["mutect"]
    shell:
        op.as_one_line("""
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
            FilterMutectCalls
            {params.opts}
            -V {input.vcf}
            -R {input.fasta}
            --tumor-segmentation {input.segments}
            --contamination-table {input.contamination}
            --ob-priors {input.model}
            --create-output-variant-index
            -O {output.vcf}
            > {log} 2>&1
        """)

# Filters for PASS and germline variants
# This will only take somatic ones if filtering just for PASSED (need to still maintain germline ones)
rule _purecn_mutect2_filter_vcf:
    input:
        vcf = str(rules._purecn_annotate_vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_passed.vcf.gz",
        tbi = CFG["dirs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_passed.vcf.gz.tbi"
    params:
        filter_for_opts = CFG["options"]["mutect2"]["filter_for"],
        filter_out_opts = CFG["options"]["mutect2"]["filter_out"]
    log:
        CFG["logs"]["mutect2"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_filter_vcf.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["post_vcf"]
    threads: 1
    shell:
        op.as_one_line("""
        bcftools view {params.filter_for_opts} -e "{params.filter_out_opts}" {input.vcf} -Oz -o {output.vcf} > {log} 2>&1 &&
         tabix -p vcf {output.vcf} >> {log} 2>&1
        """)

# PureCN - by extension rsamtools, does not have CRAM compatibility, even with R v4.5.3
# https://github.com/Bioconductor/Rsamtools/issues/21
# https://github.com/Bioconductor/Rsamtools/issues/56
# Work around is to use GATK to also calculate coverage and then feed it into pureCN for GC normalization

rule _purecn_gatk_depthOfCoverage:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        bai = str(rules._purecn_input_bam.output.bai),
        crai = str(rules._purecn_input_bam.output.crai),
        intervals = str(rules._purecn_gatk_interval_list_chrom.output.chrom_int),
        fasta = str(rules._purecn_symlink_fasta.output.fasta)
    output:
        coverage = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.sample_interval_summary",
        statistics = temp(CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}.sample_interval_statistics")
    resources:
        **CFG["resources"]["gatk_depth"]
    threads: 1
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8),
        opts = CFG["options"]["coverage"]["depth_coverage"],
        base_name = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.{chrom}",
        temp_dir = config["lcr-modules"]["_shared"]["temp_directory"]
    log:
        CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/gatk_coverage_{chrom}.log"
    conda:
        CFG["conda_envs"]["mutect"]
    shell:
        op.as_one_line(
        """
        gatk --java-options "-Xmx{params.mem_mb}m -Djava.io.tmpdir={params.temp_dir}"
            DepthOfCoverage
            {params.opts}
            --omit-depth-output-at-each-base
            --omit-locus-table
            --omit-per-sample-statistics
            --interval-merging-rule OVERLAPPING_ONLY
            -R {input.fasta}
            -I {input.bam}
            -O {params.base_name}
            -L {input.intervals}
            > {log} 2>&1
        """
        )


def _get_gatk_chr_cov_depth(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    coverage = expand(
        str(rules._purecn_gatk_depthOfCoverage.output.coverage),
        chrom = chrs,
        allow_missing = True
    )
    return(coverage)

def _get_gatk_chr_cov_stats(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    with open(checkpoints._purecn_input_chroms_withY.get(**wildcards).output.txt) as f:
        chrs = f.read().rstrip("\n").split("\n")
    statistics = expand(
        str(rules._purecn_gatk_depthOfCoverage.output.statistics),
        chrom = chrs,
        allow_missing = True
    )
    return(statistics)

rule _purecn_gatk_coverage_concatenate_depths:
    input:
        depth = _get_gatk_chr_cov_depth,
        statistics = _get_gatk_chr_cov_stats,
    output:
        depth = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.sample_interval_summary.gz"
    threads: 1
    shell:
        op.as_one_line(
            """
            file1=$(echo {input.depth} | cut -d " " -f1 ) ;
            head -n 1 $file1 | gzip > {output.depth} ;
            for sample in {input.depth};
            do
                awk '(NR > 1)' $sample | gzip >> {output.depth} ;
            done
            """)


# pureCN Coverage.R used to normalize by GC, since depth is calculated above by GATK
# These pureCN coverage files are used in making the PureCN DB that is used in the denovo case
rule _purecn_coverage:
    input:
        bam = str(rules._purecn_input_bam.output.bam),
        bai = str(rules._purecn_input_bam.output.bai),
        crai = str(rules._purecn_input_bam.output.crai),
        intervals = str(rules._purecn_symlink_intervals.output.intervals),
        coverage = str(rules._purecn_gatk_coverage_concatenate_depths.output.depth),
        fasta = str(rules._purecn_symlink_fasta.output.fasta)
    output:
        coverage = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_coverage_loess.txt.gz"
    params:
        coverage_script = CFG["software"]["coverage_script"],
        outdir = CFG["dirs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}",
        force = CFG["options"]["coverage"]["force"],
        opt = CFG["options"]["coverage"]["opts"],
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["coverage"]
    threads:
        CFG["threads"]["coverage"]
    log:
        CFG["logs"]["coverage"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}_coverage_loess.log"
    shell:
        op.as_one_line("""
            PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/ ;
            export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
            echo -e "Using {params.coverage_script} instead of default $PURECN/Coverage.R..."> {log} 2>&1;
            Rscript --vanilla {params.coverage_script}  --out-dir {params.outdir}
            --bam {input.bam}
            --name {wildcards.tumour_id}
            --reference {input.fasta}
            --coverage {input.coverage}
            --intervals {input.intervals} {params.force} {params.opt} >> {log} 2>&1
        """)


##### Run pureCN using CNVkit mode
# It will first try to use segmentation_method "Hclust" as the tool recommends
# If that fails, it will try "PSCBS"
# If that fails, it will use "none"
# Need to make the directories first so that they exist when the Rscript tries to run
rule _purecn_make_outdirs_cnvkit_mode:
    input:
        cnr = str(rules._purecn_symlink_cnvkit_cnr.output.cnr),
        seg = str(rules._purecn_symlink_cnvkit_seg.output.seg),
        vcf = str(rules._purecn_mutect2_filter_vcf.output.vcf),
        mapping_bias = str(rules._purecn_symlink_database_cnvkit.output.mapping_bias),
        blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist)
    output:
        done_pscbs = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/make_dir.done",
        done_Hclust = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/make_dir.done",
        done_none = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/make_dir.done"
    shell:
        op.as_one_line("""
        touch {output.done_pscbs};
        touch {output.done_Hclust};
        touch {output.done_none}
        """)

checkpoint _purecn_cnvkit_mode_run:
    input:
        cnr = str(rules._purecn_symlink_cnvkit_cnr.output.cnr),
        seg = str(rules._purecn_symlink_cnvkit_seg.output.seg),
        vcf = str(rules._purecn_mutect2_filter_vcf.output.vcf),
        mapping_bias = str(rules._purecn_symlink_database_cnvkit.output.mapping_bias),
        blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist),
        done_pscbs = str(rules._purecn_make_outdirs_cnvkit_mode.output.done_pscbs),
        done_Hclust = str(rules._purecn_make_outdirs_cnvkit_mode.output.done_Hclust),
        done_none = str(rules._purecn_make_outdirs_cnvkit_mode.output.done_none)
    output:
        done = CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/done"
    params:
        genome = lambda w: config["lcr-modules"]["purecn"]["options"]["genome_builds_map"][w.genome_build],
        min_offtarget = CFG["options"]["cnvkit_mode"]["min_offtarget"],
        max_cn = CFG["options"]["cnvkit_mode"]["max_cn"],
        model = CFG["options"]["cnvkit_mode"]["model"],
        alpha = CFG["options"]["cnvkit_mode"]["alpha"],
        opts = CFG["options"]["cnvkit_mode"]["opts"]
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["cnvkit_mode"]
    threads:
        CFG["threads"]["cnvkit_mode"]
    log:
        CFG["logs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        """
        set +e
        echo $CONDA_DEFAULT_ENV > {log};
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/
        export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
        echo "Trying Hclust" >> {log};
        Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_Hclust}) \
            --sampleid {wildcards.tumour_id} \
            --tumor {input.cnr} \
            --seg-file {input.seg} \
            --vcf {input.vcf} \
            --snp-blacklist {input.blacklist} \
            --mapping-bias-file {input.mapping_bias} \
            --genome {params.genome} \
            --fun-segmentation Hclust \
            --max-copy-number {params.max_cn} \
            --min-fraction-offtarget {params.min_offtarget} \
            --alpha {params.alpha} \
            --model {params.model} \
            --cores {threads} \
            {params.opts}  &>> {log}
        if [[ $? -ne 0 ]]; then
            echo "Trying PSCBS" >> {log};
            Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_pscbs}) \
                --sampleid {wildcards.tumour_id} \
                --tumor {input.cnr} \
                --seg-file {input.seg} \
                --vcf {input.vcf} \
                --snp-blacklist {input.blacklist} \
                --mapping-bias-file {input.mapping_bias} \
                --genome {params.genome} \
                --fun-segmentation PSCBS \
                --max-copy-number {params.max_cn} \
                --min-fraction-offtarget {params.min_offtarget} \
                --alpha {params.alpha} \
                --model {params.model} \
                --cores {threads} \
                {params.opts} &>> {log}
            if [[ $? -ne 0 ]]; then
                echo "Trying segmentation_method: none" >> {log};
                Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_none}) \
                    --sampleid {wildcards.tumour_id} \
                    --tumor {input.cnr} \
                    --seg-file {input.seg} \
                    --vcf {input.vcf} \
                    --snp-blacklist {input.blacklist} \
                    --mapping-bias-file {input.mapping_bias} \
                    --genome {params.genome} \
                    --fun-segmentation none \
                    --max-copy-number {params.max_cn} \
                    --min-fraction-offtarget {params.min_offtarget} \
                    --alpha {params.alpha} \
                    --model {params.model} \
                    --cores {threads} \
                    {params.opts}  &>> {log}
                if [[ $? -ne 0 ]]; then
                    exit 1
                else
                    touch {output.done}
                fi
            else
                touch {output.done}
            fi
        else
            touch {output.done}
        fi
        """

def _get_cnvkit_mode_run_result(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    base_dir = os.path.dirname(str(checkpoints._purecn_cnvkit_mode_run.get(**wildcards).output[0]))
    SEG_METHODS = glob_wildcards(base_dir + "/{seg_method}/{tumour_id}_dnacopy.seg").seg_method
    SEG_METHODS = [SEG_METHODS] if isinstance(SEG_METHODS, str) else SEG_METHODS # makes it a list when only one value is returned by the glob
    if any("Hclust" == s for s in SEG_METHODS):
        files = {
        "ploidy": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}_dnacopy.seg",
        "gene_cn": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}_genes.csv",
        "loh": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/Hclust/{tumour_id}.log"
        }
    elif any("PSCBS" == s for s in SEG_METHODS):
        files = {
        "ploidy": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}_dnacopy.seg",
        "gene_cn": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}_genes.csv",
        "loh": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.log"
        }
    else:
        files = {
        "ploidy": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}_dnacopy.seg",
        "gene_cn": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}_genes.csv",
        "loh": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_cnvkit"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.log"
        }
    return(files)

rule _purecn_symlink_cnvkit_mode_run_result:
    input:
        unpack(_get_cnvkit_mode_run_result)
    output:
        ploidy = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
        seg = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        gene_cn = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_genes.csv",
        loh = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
        log = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.log"
    run:
        op.relative_symlink(input.ploidy, output.ploidy, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)
        op.relative_symlink(input.gene_cn, output.gene_cn, in_module = True)
        op.relative_symlink(input.loh, output.loh, in_module = True)
        op.relative_symlink(input.log, output.log, in_module = True)


##### Run pureCN using denovo mode, which uses it's own coverage files
# It will first try to use segmentation_method "PSCBS"
# If that fails, it will try "CBS"
# If that fails, it will use "none"
# Need to make the directories first so that they exist when the Rscript tries to run

rule _purecn_make_outdirs_denovo_mode:
    input:
        cnr = str(rules._purecn_coverage.output.coverage),
        vcf = str(rules._purecn_mutect2_filter_vcf.output.vcf),
        mapping_bias = str(rules._purecn_symlink_database_denovo.output.mapping_bias),
        normal_db = str(rules._purecn_symlink_database_denovo.output.normal_db),
        blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist),
        intervals = str(rules._purecn_symlink_intervals.output.intervals)
    output:
        done_pscbs = CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/make_dir.done",
        done_cbs = CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/make_dir.done",
        done_none = CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/make_dir.done"
    shell:
        op.as_one_line("""
        touch {output.done_pscbs};
        touch {output.done_cbs};
        touch {output.done_none}
        """)

checkpoint _purecn_denovo_mode_run:
    input:
        cnr = str(rules._purecn_coverage.output.coverage),
        vcf = str(rules._purecn_mutect2_filter_vcf.output.vcf),
        mapping_bias = str(rules._purecn_symlink_database_denovo.output.mapping_bias),
        normal_db = str(rules._purecn_symlink_database_denovo.output.normal_db),
        blacklist = str(rules._purecn_setup_blacklist_bed.output.blacklist),
        intervals = str(rules._purecn_symlink_intervals.output.intervals),
        done_pscbs = str(rules._purecn_make_outdirs_denovo_mode.output.done_pscbs),
        done_cbs = str(rules._purecn_make_outdirs_denovo_mode.output.done_cbs),
        done_none = str(rules._purecn_make_outdirs_denovo_mode.output.done_none)
    output:
        done = CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/done"
    params:
        genome = lambda w: config["lcr-modules"]["purecn"]["options"]["genome_builds_map"][w.genome_build],
        min_offtarget = CFG["options"]["denovo_mode"]["min_offtarget"],
        max_cn = CFG["options"]["denovo_mode"]["max_cn"],
        model = CFG["options"]["denovo_mode"]["model"],
        alpha = CFG["options"]["denovo_mode"]["alpha"],
        opts = CFG["options"]["denovo_mode"]["opts"]
    conda:
        CFG["conda_envs"]["purecn"]
    resources:
        **CFG["resources"]["denovo_mode"]
    threads:
        CFG["threads"]["denovo_mode"]
    log:
        CFG["logs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}.log"
    shell:
        """
        set +e
        echo $CONDA_DEFAULT_ENV > {log};
        PURECN=$CONDA_DEFAULT_ENV/lib/R/library/PureCN/extdata/
        export R_LIBS=$CONDA_DEFAULT_ENV/lib/R/library/ ;
        echo "Trying PBCBS" >> {log};
        Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_pscbs})  \
            --sampleid {wildcards.tumour_id} \
            --tumor {input.cnr} \
            --vcf {input.vcf} \
            --snp-blacklist {input.blacklist} \
            --mapping-bias-file {input.mapping_bias} \
            --normaldb {input.normal_db} \
            --genome {params.genome} \
            --intervals {input.intervals} \
            --fun-segmentation PSCBS \
            --max-copy-number {params.max_cn} \
            --min-fraction-offtarget {params.min_offtarget} \
            --alpha {params.alpha} \
            --model {params.model} \
            --cores {threads} \
            {params.opts} &>> {log}
        if [[ $? -ne 0 ]]; then
            echo "Trying CBS" >> {log};
            Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_cbs})  \
                --sampleid {wildcards.tumour_id} \
                --tumor {input.cnr} \
                --vcf {input.vcf} \
                --snp-blacklist {input.blacklist} \
                --mapping-bias-file {input.mapping_bias} \
                --normaldb {input.normal_db} \
                --genome {params.genome} \
                --intervals {input.intervals} \
                --fun-segmentation CBS \
                --max-copy-number {params.max_cn} \
                --min-fraction-offtarget {params.min_offtarget} \
                --alpha {params.alpha} \
                --model {params.model} \
                --cores {threads} \
                {params.opts} &>> {log}
            if [[ $? -ne 0 ]]; then
                echo "Trying segmentation_method: none" >> {log};
                Rscript --vanilla $PURECN/PureCN.R --out $(dirname {input.done_none})  \
                    --sampleid {wildcards.tumour_id} \
                    --tumor {input.cnr} \
                    --vcf {input.vcf} \
                    --snp-blacklist {input.blacklist} \
                    --mappingbiasfile {input.mapping_bias} \
                    --normaldb {input.normal_db} \
                    --genome {params.genome} \
                    --intervals {input.intervals} \
                    --fun-segmentation none \
                    --max-copy-number {params.max_cn} \
                    --min-fraction-offtarget {params.min_offtarget} \
                    --alpha {params.alpha} \
                    --model {params.model} \
                    --cores {threads} \
                    {params.opts} &>> {log}
                if [[ $? -ne 0 ]]; then
                    exit 1
                else
                    touch {output.done}
                fi
            else
                touch {output.done}
            fi
        else
            touch {output.done}
        fi
        """

def _get_denovo_mode_run_result(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    base_dir = os.path.dirname(str(checkpoints._purecn_denovo_mode_run.get(**wildcards).output[0]))
    SEG_METHODS = glob_wildcards(base_dir + "/{seg_method}/{tumour_id}_dnacopy.seg").seg_method
    SEG_METHODS = [SEG_METHODS] if isinstance(SEG_METHODS, str) else SEG_METHODS # makes it a list when only one value is returned by the glob
    if any("PSCBS" == s for s in SEG_METHODS):
        files = {
        "ploidy": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}_dnacopy.seg",
        "loh": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/PSCBS/{tumour_id}.log"
        }
    elif any("CBS" == s for s in SEG_METHODS): # has potential to match "PSCBS", but that would mean it passes the if statement above, so don't need to worry about that
        files = {
        "ploidy": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}_dnacopy.seg",
        "loh": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/CBS/{tumour_id}.log"
        }
    else:
        files = {
        "ploidy": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.csv",
        "seg": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}_dnacopy.seg",
        "loh": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}_loh.csv",
        "rds": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.rds"),
        "pdf": temp(CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.pdf"),
        "log": CFG["dirs"]["pureCN_denovo"] + "{seq_type}--{genome_build}/{capture_space}/{tumour_id}/none/{tumour_id}.log"
        }
    return(files)

rule _purecn_symlink_denovo_mode_run_result:
    input:
        unpack(_get_denovo_mode_run_result)
    output:
        ploidy = CFG["dirs"]["pureCN_final_result"] + "purecn_denovo/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv",
        seg = CFG["dirs"]["pureCN_final_result"] + "purecn_denovo/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        loh = CFG["dirs"]["pureCN_final_result"] + "purecn_denovo/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv",
        log = CFG["dirs"]["pureCN_final_result"] + "purecn_denovo/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.log"
    run:
        op.relative_symlink(input.ploidy, output.ploidy, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)
        op.relative_symlink(input.loh, output.loh, in_module = True)
        op.relative_symlink(input.log, output.log, in_module = True)
        
# Record segmentation method info for the final_result seg
rule _purecn_record_seg_fn:
    input:
        seg = CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
    output:
        txt = CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_segmentation_fn.txt"
    log:
        log = CFG["logs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_segmentation_fn.log"
    shell:
        op.as_one_line("""
        echo "segmentation_fn" > {output.txt};
        SEG=$(readlink -f {input.seg}) > {log.log} 2>&1;
        IFS="/" read -ra path_parts <<< $SEG >> {log.log} 2>&1;
        echo "${{path_parts[${{#path_parts[@]}}-2]}}" >> {output.txt} 2>> {log.log} 
        """)

rule _purecn_cnv2igv:
    input:
        seg = CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg",
        txt = str(rules._purecn_record_seg_fn.output.txt)
    output:
        seg = CFG["dirs"]["cnv2igv"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg"
    log:
        stderr = CFG["logs"]["cnv2igv"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/cnv2igv_{tumour_id}.stderr.log"
    params:
        cnv2igv_script = CFG["options"]["cnv2igv_script_path"],
        opts = CFG["options"]["preserve"]
    conda:
        CFG["conda_envs"]["cnv2igv"]
    threads: 1
    shell:
        op.as_one_line("""
        python {params.cnv2igv_script} --mode purecn {params.opts} --sample {wildcards.tumour_id} {input.seg} > {output.seg} 2>> {log.stderr}
        """)
        
def _purecn_get_chain(wildcards):
    if "38" in str({wildcards.genome_build}):
        return reference_files("genomes/{genome_build}/chains/grch38/hg38ToHg19.over.chain")
    else:
        return reference_files("genomes/{genome_build}/chains/grch37/hg19ToHg38.over.chain")
        

# Convert the coordinates of seg file to a different genome build
rule _purecn_convert_coordinates:
    input:
        purecn_native = str(rules._purecn_cnv2igv.output.seg),
        purecn_chain = _purecn_get_chain
    output:
        purecn_lifted = CFG["dirs"]["convert_coordinates"] + "{mode}/from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}.lifted_{chain}.seg"
    log:
        stderr = CFG["logs"]["convert_coordinates"] + "{mode}/from--{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.lifted_{chain}.stderr.log"
    threads: 1
    params:
        liftover_script = CFG["options"]["liftover_script_path"],
        liftover_minmatch = CFG["options"]["liftover_minMatch"]
    conda:
        CFG["conda_envs"]["liftover"]
    shell:
        op.as_one_line("""
        bash {params.liftover_script}
        SEG
        {input.purecn_native}
        {output.purecn_lifted}
        {input.purecn_chain}
        YES
        {params.liftover_minmatch}
        2>> {log.stderr}
        """)
        
def _purecn_prepare_projection(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_genome_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    if "38" in this_genome_build[0]:
        hg38_projection = str(rules._purecn_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_projection = str(rules._purecn_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        grch37_projection = grch37_projection.replace("{chain}", "hg38ToHg19")
    else:
        grch37_projection = str(rules._purecn_cnv2igv.output.seg).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_projection = str(rules._purecn_convert_coordinates.output.purecn_lifted).replace("{genome_build}", this_genome_build[0]).replace("{capture_space}", this_space[0])
        hg38_projection = hg38_projection.replace("{chain}", "hg19ToHg38")
    return{
        "grch37_projection": grch37_projection,
        "hg38_projection": hg38_projection
    }
    

# Fill the missing segments of seg files with neutral regions to complete the genome coverage
rule _purecn_fill_segments:
    input:
        unpack(_purecn_prepare_projection)
    output:
        grch37_filled = temp(CFG["dirs"]["fill_regions"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}.grch37.seg"),
        hg38_filled = temp(CFG["dirs"]["fill_regions"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}.hg38.seg")
    log:
        stderr = CFG["logs"]["fill_regions"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}_fill_segments.stderr.log"
    threads: 1
    params:
        path = config["lcr-modules"]["_shared"]["lcr-scripts"] + "fill_segments/" + CFG["options"]["fill_segments_version"]
    conda:
        CFG["conda_envs"]["bedtools"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
        echo "Filling grch37 projection" >> {log.stderr};
        bash {params.path}fill_segments.sh
        {params.path}src/chromArm.grch37.bed
        {input.grch37_projection}
        {params.path}src/blacklisted.grch37.bed
        {output.grch37_filled}
        {wildcards.tumour_id}
        SEG
        2>> {log.stderr};
        echo "Filling hg38 projection" >> {log.stderr};
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
rule _purecn_normalize_projection:
    input:
        filled = CFG["dirs"]["fill_regions"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}.{projection}.seg",
        chrom_file = reference_files("genomes/{projection}/genome_fasta/main_chromosomes.txt")
    output:
        projection = CFG["dirs"]["normalize"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}.{projection}.seg"
    resources:
        **CFG["resources"]["post_purecn"]
    threads: 1
    run:
        # read the main chromosomes file of the projection
        chromosomes = pd.read_csv(input.chrom_file, sep = "\t", names=["chromosome"], header=None)
        # handle chr prefix
        if "chr" in chromosomes["chromosome"][0]:
            seg_open = pd.read_csv(input.filled, sep = "\t")
            chrom = list(seg_open['chrom'])
            # avoid cases of chrchr1 if the prefix already there
            for i in range(len(chrom)):
                if 'chr' not in str(chrom[i]):
                    chrom[i]='chr'+str(chrom[i])
            seg_open.loc[:, 'chrom']=chrom
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')
        else:
            # remove chr prefix
            seg_open = pd.read_csv(input.filled, sep = "\t")
            seg_open["chrom"] = seg_open["chrom"].astype(str).str.replace('chr', '')
            seg_open.to_csv(output.projection, sep="\t", index=False, na_rep='NA')

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _purecn_output_projection:
    input:
        projection = str(rules._purecn_normalize_projection.output.projection)
    output:
        projection = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--projection/{tumour_id}.{tool}.{projection}.seg"
    wildcard_constraints:
        mode = "|".join(CFG["modes"])
    run:
        op.relative_symlink(input.projection, output.projection, in_module = True)

def _get_capture_space(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()

    inputs = {
        "raw_seg": CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_dnacopy.seg".replace("{capture_space}", this_space[0]),
        "ploidy": CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.csv".replace("{capture_space}", this_space[0]),
        "loh": CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_loh.csv".replace("{capture_space}", this_space[0]),
        "log": CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}.log".replace("{capture_space}", this_space[0]),
        "txt": CFG["dirs"]["pureCN_final_result"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}/{tumour_id}_segmentation_fn.txt".replace("{capture_space}", this_space[0]),
        "seg": CFG["dirs"]["cnv2igv"] + "{mode}/{seq_type}--{genome_build}/{capture_space}/{tumour_id}.seg".replace("{capture_space}", this_space[0])
    }
    return inputs
    
rule _purecn_output:
    input:
        unpack(_get_capture_space)
    output:
        raw_seg = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}_dnacopy.seg",
        ploidy = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}.csv",
        loh = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}_loh.csv",
        log = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}.log",
        txt = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}_segmentation_fn.txt",
        seg = CFG["dirs"]["outputs"] + "{mode}/{seq_type}--{genome_build}/{tumour_id}/{tumour_id}_cnv2igv.seg"
    run:
        op.relative_symlink(input.raw_seg, output.raw_seg, in_module = True)
        op.relative_symlink(input.ploidy, output.ploidy, in_module = True)
        op.relative_symlink(input.loh, output.loh, in_module = True)
        op.relative_symlink(input.log, output.log, in_module = True)
        op.relative_symlink(input.txt, output.txt, in_module = True)
        op.relative_symlink(input.seg, output.seg, in_module = True)

# Decide best seg out of the two modes
# Select the seg file that had the least amount of deviation (noise) - as measured by the MAD value
# If it's the de novo file, check that it's not over segmented
def _purecn_take_lowest_MAD(wildcards):
    CFG = config["lcr-modules"]["purecn"]
    tbl = CFG["runs"]
    this_build = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_genome_build"].tolist()[0]
    this_space = tbl[(tbl.tumour_sample_id == wildcards.tumour_id) & (tbl.tumour_seq_type == wildcards.seq_type)]["tumour_capture_space"].tolist()[0]
    
    cnvkit_seg = CFG["dirs"]["outputs"] + "purecn_cnvkit/" + wildcards.seq_type + "--projection/" + wildcards.tumour_id + "." + wildcards.tool + "." + wildcards.projection + ".seg"
    denovo_seg = CFG["dirs"]["outputs"] + "purecn_denovo/" + wildcards.seq_type + "--projection/" + wildcards.tumour_id + "." + wildcards.tool + "." + wildcards.projection + ".seg"

    cnvkit_log = CFG["dirs"]["pureCN_final_result"] + "purecn_cnvkit/" + wildcards.seq_type + "--" + this_build + "/" + this_space + "/" + wildcards.tumour_id + "/" + wildcards.tumour_id + ".log"
    denovo_log = CFG["dirs"]["pureCN_final_result"] + "purecn_denovo/" + wildcards.seq_type + "--" + this_build + "/" + this_space + "/" + wildcards.tumour_id + "/" + wildcards.tumour_id + ".log"

    if os.path.isfile(str(cnvkit_seg)) and os.path.isfile(str(denovo_seg)):

        cnvkit_mapd = list()
        with open(cnvkit_log) as file:
            for line in file.readlines():
                if 'Mean standard deviation of log-ratios' in line:
                    cnvkit_mapd.append(str(line))

        # Take the 8th element - representing the MAD score
        cnvkit_mapd_value = float(cnvkit_mapd[len(cnvkit_mapd)-1].split()[8])

        denovo_mapd = list()
        with open(denovo_log) as file:
            for line in file.readlines():
                if 'Mean standard deviation of log-ratios' in line:
                    denovo_mapd.append(str(line))

        # Take the 8th element - representing the MAD score
        denovo_mapd_value = float(denovo_mapd[len(denovo_mapd)-1].split()[8])

        if (denovo_mapd_value < cnvkit_mapd_value):
            # Account for seg files that are over-segmented but have a lower MAD
            count_denovo = []
            with open(denovo_seg) as file:
                next(file)
                for line in file:
                    line = line.rstrip('\n').rstrip('\r')
                    cols = line.split('\t')
                    if (float(cols[7]) > 0.5 or float(cols[7]) <-0.5 and float(cols[8]) != 1):
                        count_denovo.append(line)
            count_denovo = len(count_denovo)

            count_cnvkit = []
            with open(cnvkit_seg) as file:
                next(file)
                for line in file:
                    line = line.rstrip('\n').rstrip('\r')
                    cols = line.split('\t')
                    if (float(cols[7]) > 0.5 or float(cols[7]) <-0.5 and float(cols[8]) != 1):
                        count_cnvkit.append(line)
            count_cnvkit = len(count_cnvkit)
            if (count_cnvkit > 0):
                if (float(count_denovo/count_cnvkit) > 2.8 and count_denovo > 400):
                    best_seg = str(cnvkit_seg)
                else:
                    best_seg = str(denovo_seg)
            else:
                best_seg = str(denovo_seg)
        else:
            best_seg = str(cnvkit_seg)

    elif os.path.isfile(str(cnvkit_seg)):
        best_seg = cnvkit_seg
    else:
        best_seg = denovo_seg
        
    return{
        "best_seg": best_seg
    }
    
rule _purecn_best_seg:
    input:
        unpack(_purecn_take_lowest_MAD)
    output:
        best_seg = CFG["dirs"]["outputs"] + "best_seg/{seq_type}--projection/{tumour_id}.{tool}.{projection}.seg"
    run:
        op.relative_symlink(input.best_seg, output.best_seg, in_module = True)

rule _purecn_best_seg_info:
    input:
        seg = str(rules._purecn_best_seg.output.best_seg)
    output:
        tsv = CFG["dirs"]["outputs"] + "best_seg/{seq_type}--projection/{tumour_id}.{tool}.{projection}.info.tsv"
    log:
        log = CFG["logs"]["outputs"] + "best_seg/{seq_type}--projection/{tumour_id}.{tool}.{projection}.info.log"
    shell:
        op.as_one_line("""
        echo -e "winning_mode" > {output.tsv};
        SEG=$(readlink -f {input.seg}) > {log.log} 2>&1;
        IFS="/" read -ra path_parts <<< $SEG >> {log.log} 2>&1;
        echo -e "${{path_parts[${{#path_parts[@]}}-3]}}" >> {output.tsv} 2>> {log.log}
        """)
    
rule _purecn_all:
    input:
        expand(
            expand(
            [
                str(rules._purecn_record_seg_fn.output.txt),
                str(rules._purecn_output.output.raw_seg),
                str(rules._purecn_output.output.ploidy),
                str(rules._purecn_output.output.loh),
                str(rules._purecn_output.output.log),
                str(rules._purecn_output.output.seg)
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            capture_space=CFG["runs"]["tumour_capture_space"],
            allow_missing = True
            ),
        mode = CFG["modes"]
        ),
        expand(
            expand(
            [
                str(rules._purecn_output_projection.output.projection)
            ],
            zip, # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            allow_missing = True
            ),
        tool = "purecn",
        projection=CFG["requested_projections"],
        mode = CFG["modes"]
        ),
        expand(
            expand(
            [
                str(rules._purecn_best_seg.output.best_seg),
                str(rules._purecn_best_seg_info.output.tsv)
            ],
            zip, # Run expand() with zip(), not product()
            tumour_id=CFG["runs"]["tumour_sample_id"],
            seq_type=CFG["runs"]["tumour_seq_type"],
            allow_missing = True
            ),
        tool = "purecn",
        projection=CFG["requested_projections"]
        )
                


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
