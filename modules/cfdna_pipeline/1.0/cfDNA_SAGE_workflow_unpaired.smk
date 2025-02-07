#!/usr/bin/env snakemake

import os
import pandas as pd
import snakemake
import glob
snakemake.utils.min_version("7")

import sys
MODULE_PATH = os.path.join(config["lcr-modules"]["__shared"]["lcr-modules"], "/modules/cfdna_pipeline/1.0/")
sys.path.append(MODULE_PATH) # add local module to path

# generate paths for file locations
BAM_OUTDIR = os.path.join(config["lcr-modules"]["__shared"]["root_output_dir"], "bam_pipeline")
UTILSDIR = os.path.join(config["lcr-modules"]["__shared"]["lcr-modules"], "/modules/cfdna_pipeline/1.0/utils")
SAGE_OUTDIR = os.path.join(config["lcr-modules"]["__shared"]["root_output_dir"], "sage_pipeline")

all_samples = config["lcr-modules"]["_shared"]["samples"].copy()
SAMPLESHEET_UN = all_samples.loc[(all_samples["tissue_status"] == "tumor") & (all_samples["matched_normal"] == "unpaired")].copy()

####################################### input functions

def older_sample_mafs(wildcards):
    """ Return a list of older sample mafs for augmenting the current maf file.
    """
    # get patient_id
    patient_id = unmatched_samplesheet.loc[unmatched_samplesheet["sample_id"] == wildcards.sample, "patient_id"].values[0]
    # get all samples for this patient
    patient_samples = unmatched_samplesheet.loc[(unmatched_samplesheet["patient_id"] == patient_id) & (unmatched_samplesheet['timepoint'] != 'normal' )]["sample_id"].tolist()
    patient_samples = [s for s in patient_samples if s != wildcards.sample]

    return expand(os.path.join(SAGE_OUTDIR, "12-filtered_unmatched/{sample}.sage.filtered.unmatched.maf"), sample=patient_samples)

########################################################### Run variant calling
rule run_sage_un:
    input:
        tbam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam"),
    output:
        vcf = os.path.join(SAGE_OUTDIR, "01-SAGE/{sample}/{sample}.unmatched.sage.vcf")
    params:
        ref_genome = config["ref_genome"],
        ref_genome_version = "38" if config["ref_genome_ver"] == "GRCh38" else "37",
        # Panel regions and hotspots
        hotspots_vcf = config["cappseq_snv_pipeline"]["sage_hotspots"],
        panel_regions = config["captureregions"],
        ensembl = config['cappseq_snv_pipeline']['ensembl'],
        # Miscellaneous
        max_depth = config["cappseq_snv_pipeline"]["max_depth"],
        min_map = config["cappseq_snv_pipeline"]["min_map_qual"],
        hard_vaf_cutoff = config["cappseq_snv_pipeline"]["tumor_min_vaf"],
        soft_vaf_cutoff = config["cappseq_snv_pipeline"]["tumor_soft_min_vaf"],
        min_norm_depth = 7,
        alt_support = config["cappseq_snv_pipeline"]["min_alt_depth"]
    resources:
        mem_mb = 10000,
        runtime_min = 60
    threads: 12
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.sage_run.log")
    conda:
        "envs/sage.yaml"
    shell:
        """sage -tumor {wildcards.sample} -tumor_bam {input.tbam} \
        -output_vcf {output.vcf} -ref_genome {params.ref_genome} \
        -ref_genome_version {params.ref_genome_version} \
        -hotspots {params.hotspots_vcf} -panel_bed {params.panel_regions} \
        -high_confidence_bed {params.panel_regions} \
        -ensembl_data_dir {params.ensembl} \
        -max_read_depth {params.max_depth} -min_map_quality {params.min_map} \
        -hard_min_tumor_vaf {params.hard_vaf_cutoff} -hard_min_tumor_raw_alt_support {params.alt_support} \
        -panel_min_tumor_vaf {params.soft_vaf_cutoff} \
        -high_confidence_min_tumor_vaf {params.soft_vaf_cutoff} \
        -bqr_min_map_qual {params.min_map} -threads {threads} &> {log}
        """

rule filter_sage_un:
    input:
        vcf = rules.run_sage_un.output.vcf
    output:
        vcf = os.path.join(SAGE_OUTDIR, "02-vcfs/{sample}.sage.unmatched.passed.vcf")
    conda:
        "envs/bcftools.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    group: "filter_sage"
    params:
        notlist = config["cappseq_snv_pipeline"]["notlist"]
    shell:
        """
        bcftools view -T ^{params.notlist} -f PASS {input.vcf} -O vcf -o {output.vcf}
        """

# Flag positions with a high incidence of masked bases
rule flag_masked_pos_un:
    input:
        bam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam")
    output:
        bed_raw = temp(os.path.join(SAGE_OUTDIR, "03-masked_pos/{sample}.maskedpos.bed")),
        bed = os.path.join(SAGE_OUTDIR, "03-masked_pos/{sample}.maskedpos.bed.gz")
    params:
        script = os.path.join(UTILSDIR, "mask_n_sites.py"),
        n_threshold = config["cappseq_snv_pipeline"]["mask_threshold"],
        min_count = config["cappseq_snv_pipeline"]["mask_count"],
        panel_regions = config["captureregions"]
    conda:
        "envs/bcftools.yaml"
    resources:
        mem_mb = 10000
    threads: 1
    group: "filter_sage"
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.maskpos.log")
    shell:
        """
        {params.script} --input {input.bam} --regions {params.panel_regions} --output {output.bed_raw} --count {params.min_count} --fraction {params.n_threshold} > {log} &&
        bgzip -c {output.bed_raw} > {output.bed} && tabix -p bed {output.bed} >> {log}
        """

rule restrict_to_capture_un:
    input:
        vcf = rules.filter_sage_un.output.vcf,
        bed = rules.flag_masked_pos_un.output.bed
    output:
        vcf = os.path.join(SAGE_OUTDIR, "04-capturespace/{sample}.capspace.unmatched.vcf")
    params:
        panel_regions = config["captureregions"]
    conda:
        "envs/bcftools.yaml"
    resources:
        mem_mb = 10000
    threads: 2
    group: "filter_sage"
    shell:
        """
        bedtools intersect -a {input.vcf} -header -b {params.panel_regions} | bedtools intersect -a - -header -b {input.bed} -v | awk -F '\\t' '$4 !~ /N/ && $5 !~ /N/' | bcftools norm -m +any -O vcf -o {output.vcf}
        """

# Generate review BAM files containing only the reads supporting these variants
rule review_consensus_reads_un:
    input:
        bam_cons = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam"),
        bam_uncons = os.path.join(BAM_OUTDIR, "04-umigrouped", "{sample}.umigrouped.sort.bam"),
        vcf = rules.restrict_to_capture_un.output.vcf
    output:
        tmp_sort = temp(os.path.join(SAGE_OUTDIR, "06-supportingreads/{sample}/{sample}.umigrouped.sort.um.bam")),
        tmp_index = temp(os.path.join(SAGE_OUTDIR, "06-supportingreads/{sample}/{sample}.umigrouped.sort.bam.um.bai")),
        consensus_bam = os.path.join(SAGE_OUTDIR, "06-supportingreads/{sample}/{sample}.consensus.um.bam"),
        grouped_bam = os.path.join(SAGE_OUTDIR, "06-supportingreads/{sample}/{sample}.grouped.um.bam")
    resources:
        mem_mb = 10000
    threads: 2
    group: "filter_sage"
    params:
        ref_genome = config["ref_genome"],
        outdir = os.path.join(SAGE_OUTDIR, "06-supportingreads/{sample}/")
    conda:
        "envs/fgbio.yaml"
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.reviewconsensusvariant.log")
    shell:
        """
        samtools sort -@ 2 {input.bam_uncons} > {output.tmp_sort} && samtools index -@ 2 {output.tmp_sort} &&
        fgbio ReviewConsensusVariants --input {input.vcf} --consensus {input.bam_cons} --grouped-bam {output.tmp_sort} --ref {params.ref_genome} --output {params.outdir}/{wildcards.sample} --sample {wildcards.sample} 2>&1 > {log}
        """

rule vcf2maf_annotate_un:
    input:
        vcf = rules.restrict_to_capture_un.output.vcf
    output:
        vep_vcf = temp(os.path.join(SAGE_OUTDIR, "04-capturespace/{sample}.capspace.unmatched.vep.vcf")),
        maf = os.path.join(SAGE_OUTDIR, "05-MAFs/{sample}.sage.unmatched.maf")
    params:
        custom_enst = os.path.join(config["repo_path"], "resources/custom_enst.hg38.txt"),
        vep_data = config["cappseq_snv_pipeline"]["vep_data"],
        centre = config["cappseq_snv_pipeline"]["centre"],
        ref_fasta = config["ref_genome"],
        ref_ver = config["ref_genome_ver"]
    group: "filter_sage"
    resources:
        mem_mb = 5000
    threads: 2
    conda:
        "envs/vcf2maf.yaml"
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.vcf2maf.log")
    shell:
        """
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} \
        --tumor-id {wildcards.sample} \
        --vep-path $CONDA_PREFIX/bin/ \
        --vep-data {params.vep_data} --vep-forks {threads} \
        --custom-enst {params.custom_enst} --ref-fasta {params.ref_fasta} \
        --species homo_sapiens --ncbi-build {params.ref_ver} \
        --retain-info LPS,LPS_RC \
        --maf-center {params.centre} 2> {log} >> {log}
        """

rule augment_ssm_ChrisIndels_un:
    input:
        maf = rules.vcf2maf_annotate_un.output.maf,
        tbam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam")
    output:
        maf = os.path.join(SAGE_OUTDIR, "08-augmentssm/{sample}.sage.augment.unmatched.maf")
    resources:
        mem_mb = 10000
    threads: 4
    conda:
        "envs/augment_ssm.yaml"
    group: "filter_sage"
    params:
        script = os.path.join(UTILSDIR, "augment_ssm_ChrisIndels.py"),
        ref_genome_version = config["ref_genome_ver"]
    shell: """
    python {params.script} \
    --bam {input.tbam} \
    --maf {input.maf} \
    --genome {params.ref_genome_version} \
    --output {output.maf} \
    --threads {threads}
    """

rule filter_repetitive_seq_un:
    """
    Remove mutations which are adjacent to a repetitive sequence.
    """
    input:
        maf = rules.augment_ssm_ChrisIndels_un.output.maf
    output:
        maf = os.path.join(SAGE_OUTDIR, "10-filter_repeat/{sample}.sage.repeat_filt.umatched.maf")
    params:
        max_repeat_len = 6,
        ref_fasta = config["ref_genome"]
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    group: "filter_sage"
    shell:
        f"""python {os.path.join(UTILSDIR, "filter_rep_seq.py")} \
        --in_maf {{input.maf}} --out_maf {{output.maf}} --reference_genome {{params.ref_fasta}} --max_repeat_len {{params.max_repeat_len}}
        """

rule add_UMI_support_un:
    input:
        maf = rules.filter_repetitive_seq_un.output.maf,
        umi_bam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam")
    output:
        omaf = os.path.join(SAGE_OUTDIR, "11-umi_unmatched/{sample}.umi_support.unmatched.maf")
    params:
        script = os.path.join(UTILSDIR, "FetchVariantUMIs.py"),
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    group: "filter_sage"
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.add_UMI_support.log")
    shell:
        f"""python {{params.script}} --input_maf {{input.maf}} --input_bam {{input.umi_bam}} --output_maf {{output.omaf}} &> {{log}}
        """

rule custom_filters_un:
    input:
        maf = rules.add_UMI_support_un.output.omaf
    output:
        maf = temp(os.path.join(SAGE_OUTDIR, "12-filtered_unmatched/{sample}.sage.filtered.unmatched.maf"))
    params:
        exac_freq = float(config["cappseq_snv_pipeline"]["exac_max_freq"]),
        script = os.path.join(UTILSDIR, "unmatched_filters.py"),
        hotspot_txt = config["hotspot_manifest"],
        blacklist_txt = config["blacklist_manifest"],
        min_alt_depth = config["cappseq_snv_pipeline"]["min_alt_depth"]
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.custom_filters.log")
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    group: "filter_sage"
    shell:
        f"""python {{params.script}} --input_maf {{input.maf}} --output_maf {{output.maf}} \
        --min_alt_depth_tum {{params.min_alt_depth}} \
        --blacklist {{params.blacklist_txt}} --hotspots {{params.hotspot_txt}} --gnomad_threshold {{params.exac_freq}} &> {{log}}
        """

rule augment_maf_un:
    input:
        index_maf = rules.custom_filters_un.output.maf,
        additional_mafs = older_sample_mafs,
        index_bam = os.path.join(BAM_OUTDIR, "99-final", "{sample}.consensus.mapped.annot.bam")
    output:
        maf = os.path.join(SAGE_OUTDIR, "99-final_unmatched/{sample}.processed.unmatched.maf")
    params:
        script = os.path.join(UTILSDIR, "augmentMAF.py"),
        ref_genome_version = config["ref_genome_ver"],
        alt_support = config["cappseq_snv_pipeline"]["min_alt_depth"]
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.augment_maf.log")
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 20000
    threads: 3
    shell:
        f"""python {{params.script}} --sample_id {{wildcards.sample}} --index_maf {{input.index_maf}} --index_bam {{input.index_bam}} \
            --add_maf_files {{input.additional_mafs}} --genome_build {{params.ref_genome_version}} --threads {{threads}} \
        --alt_count_min {{params.alt_support}} --output {{output.maf}} &> {{log}}
        """

rule igv_screenshot_variants_un:
    input:
        maf = rules.augment_maf_un.output.maf,
        tbam = os.path.join(BAM_OUTDIR, "99-final/", "{sample}.consensus.mapped.annot.bam"),
    output:
        html = os.path.join(SAGE_OUTDIR, "07-IGV/{sample}_report.unmatched.html")
    params:
        refgenome = config["ref_genome"],
        cytoband = config["cappseq_snv_pipeline"]["cytoband"],
        genes = config["cappseq_snv_pipeline"]["genetrack"],
        sample_name = lambda w: w.sample
    conda:
        "envs/igv.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.igv.log")
    shell:
        """
        create_report --fasta {params.refgenome} --type mutation --tracks {input.tbam} {params.genes} --flanking 1500 --output {output.html} --standalone --title {params.sample_name} --ideogram {params.cytoband} {input.maf} > {log}
        """

rule all_sage:
    input:
        expand(str(rules.augment_maf_un.output.maf), sample=SAMPLESHEET_UN.loc[SAMPLESHEET_UN["tissue_status"]=="tumor"]["sample_id"]),
        expand(str(rules.igv_screenshot_variants_un.output.html), sample = SAMPLESHEET_UN.loc[SAMPLESHEET_UN["tissue_status"]=="tumor"]["sample_id"])
