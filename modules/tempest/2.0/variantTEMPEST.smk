#!/usr/bin/env snakemake
import os
import pandas as pd
import snakemake
import glob
snakemake.utils.min_version("7")

import sys
MODULE_PATH = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "modules/tempest/1.0/")
sys.path.append(MODULE_PATH) # add local module to path
from _version import __version__ as pv # get pipeline version

# generate paths for file locations
BAM_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "tempest/2.0/", "bam_pipeline")
UTILSDIR = os.path.join(MODULE_PATH, "utils")
SAGE_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "tempest/2.0/", "sage_pipeline")

all_samples = config["lcr-modules"]["_shared"]["samples"]
# make sure no unmatched samples are fed into workflow
SAMPLESHEET = all_samples.loc[all_samples["matched_normal"] != "unmatched"]

# compensate for canadian spelling of tumour in tissue_status
SAMPLESHEET.loc[SAMPLESHEET["tissue_status"] == "tumour", "tissue_status"] = "tumor"

##################### input functions ############################
def get_normal_bam(wildcards):
    """ Return the path to the normal bam file for the given sample.

    Input function for rules in in sage pipeline. Automatically finds
    lates version of samplesheet created by checkpoint after the bams
    are made.
    """
    normal_sample = SAMPLESHEET.loc[SAMPLESHEET["sample_id"] == wildcards.sample, "matched_normal"]

    # check if normal sample is found
    if normal_sample.empty:
        print(f"############# Normal sample for {wildcards.sample} is {normal_sample}")
        raise ValueError(f"Normal sample not found for {wildcards.sample}")
    else:
        normal_sample = normal_sample.values[0]

    # if file not found and needs to be made, return path where pipeline can make one
    return expand(config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"], sample=normal_sample)

def get_normal_name(wildcards):
    return SAMPLESHEET.loc[SAMPLESHEET["sample_id"] == wildcards.sample, "matched_normal"].values[0]

def older_sample_mafs(wildcards):
    """ Return a list of older sample mafs for augmenting the current maf file.
    """
    # get patient_id
    patient_id = SAMPLESHEET.loc[SAMPLESHEET["sample_id"] == wildcards.sample, "patient_id"].values[0]
    # get all samples for this patient
    patient_samples = SAMPLESHEET.loc[(SAMPLESHEET["patient_id"] == patient_id) & (SAMPLESHEET['tissue_status'] != 'normal' )]["sample_id"].values
    # remove the current sample
    patient_samples = [s for s in patient_samples if s != wildcards.sample]

    return expand(os.path.join(SAGE_OUTDIR, "12-filtered/{sample}.sage.filtered.maf"), sample=patient_samples)

def get_capture_space(wildcards):
    """Get the capture regions file path for a sample from the sample sheet and config"""
    capture_space = SAMPLESHEET.loc[SAMPLESHEET["sample_id"] == wildcards.sample, "capture_space"]
    
    # check if capture space is found in sample sheet
    if capture_space.empty:
        raise ValueError(f"Capture space not found for {wildcards.sample}")
    
    capture_space_value = capture_space.values[0]
    
    # check if capture space is found in config
    if capture_space_value not in config["lcr-modules"]["_shared"]["captureregions"]:
        raise ValueError(f"Capture space '{capture_space_value}' not found in config")
        
    return config["lcr-modules"]["_shared"]["captureregions"][capture_space_value]


########################################################### Run variant calling

localrules:
    filter_and_restrict_sage,
    custom_filters

def sage_dynamic_mem(wildcards, attempt, input):
    if attempt == 1:
        return max(10000, (2 * input.size_mb + 2000))
    elif attempt == 2:
        return max(15000, (input.size_mb * 3 + 5000))
    elif attempt == 3:
        return max(20000,(input.size_mb * 4))

def sage_java_mem(wildcards, attempt, input):
    if attempt == 1:
        return max(10000, (2 * input.size_mb + 2000))
    elif attempt == 2:
        return (input.size_mb * 3 + 5000) -2000
    elif attempt == 3:
        return (input.size_mb * 4) - 2000

rule run_sage:
    input:
        tbam = config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"],
        nbam = get_normal_bam
    output:
        vcf = os.path.join(SAGE_OUTDIR, "01-SAGE/{sample}/{sample}.sage.vcf")
    params:
        ref_genome = config["lcr-modules"]["_shared"]["ref_genome"],
        ref_genome_version = "38" if config["lcr-modules"]["_shared"]["ref_genome_ver"] == "GRCh38" else "37",
        # Panel regions and hotspots and inputs
        hotspots_vcf = config["lcr-modules"]["cfDNA_SAGE_workflow"]["hotspot_vcf"],
        panel_regions = lambda wildcards: get_capture_space(wildcards),
        high_conf_bed = config["lcr-modules"]["cfDNA_SAGE_workflow"]["high_conf_bed"],
        ensembl = config["lcr-modules"]["cfDNA_SAGE_workflow"]['ensembl'],
        normal_name = get_normal_name,
        # soft filters
        panel_vaf_threshold = config["lcr-modules"]["cfDNA_SAGE_workflow"]["tumor_panel_min_vaf"],
        min_norm_depth = config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_germline_depth"],
        pan_max_germ_rel_raw_bq = config["lcr-modules"]["cfDNA_SAGE_workflow"]["panel_max_germ_rel_raw_bq"],
        hotspot_max_germ_rel_raw_bq = config["lcr-modules"]["cfDNA_SAGE_workflow"]["hotspot_max_germ_rel_raw_bq"],
        hotspot_max_germline_vaf = config["lcr-modules"]["cfDNA_SAGE_workflow"].get("hotspot_max_germline_vaf", 0.1),
        # hard filters
        hard_vaf_cutoff = config["lcr-modules"]["cfDNA_SAGE_workflow"]["hard_min_vaf"],
        max_alt_norm_depth = config["lcr-modules"]["cfDNA_SAGE_workflow"]["max_normal_alt_depth"],
        min_map = config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_map_qual"],
        max_depth = config["lcr-modules"]["cfDNA_SAGE_workflow"]["max_depth"],
    resources:
        mem_mb = sage_dynamic_mem,
        java_mem = sage_java_mem,
        runtime_min = 60
    threads: 12
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.sage_run.log")
    conda:
        "envs/sage.yaml"
    shell:
        """sage -Xmx{resources.java_mem}m -tumor {wildcards.sample} -tumor_bam {input.tbam} \
        -output_vcf {output.vcf} -ref_genome {params.ref_genome} \
        -ref_genome_version {params.ref_genome_version} \
        -reference {params.normal_name} \
        -reference_bam {input.nbam} \
        -ensembl_data_dir {params.ensembl} \
        -hotspots {params.hotspots_vcf} \
        -panel_bed {params.panel_regions} \
        -high_confidence_bed {params.high_conf_bed} \
        -panel_min_tumor_vaf {params.panel_vaf_threshold} \
        -hotspot_min_tumor_vaf {params.panel_vaf_threshold} \
        -panel_min_germline_depth {params.min_norm_depth} \
        -hotspot_min_germline_depth {params.min_norm_depth} \
        -panel_max_germline_rel_qual {params.pan_max_germ_rel_raw_bq} \
        -hotspot_max_germline_rel_qual {params.hotspot_max_germ_rel_raw_bq} \
        -hotspot_max_germline_vaf {params.hotspot_max_germline_vaf} \
        -min_map_quality {params.min_map} \
        -hard_min_tumor_vaf {params.hard_vaf_cutoff} \
        -filtered_max_germline_alt_support {params.max_alt_norm_depth} \
        -max_read_depth {params.max_depth} \
        -bqr_min_map_qual {params.min_map} \
        -skip_msi_jitter \
        -threads {threads} &> {log}
        """

rule filter_and_restrict_sage:
    input:
        vcf = rules.run_sage.output.vcf,
    output:
        vcf = os.path.join(SAGE_OUTDIR, "04-capturespace/{sample}.capspace.vcf")
    params:
        panel_regions = lambda wildcards: get_capture_space(wildcards),
        notlist = config["lcr-modules"]["cfDNA_SAGE_workflow"]["notlist"]
    conda:
        "envs/bcftools.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        """
        # First apply PASS filter and optional notlist
        if [ -n "{params.notlist}" ] && [ -f "{params.notlist}" ]; then
            bcftools view -T ^{params.notlist} -f PASS {input.vcf} -O vcf | \
            bedtools intersect -a - -header -b {params.panel_regions} | \
            awk -F '\\t' '$4 !~ /N/ && $5 !~ /N/' | \
            bcftools norm -m +any -O vcf -o {output.vcf}
        else
            bcftools view -f PASS {input.vcf} -O vcf | \
            bedtools intersect -a - -header -b {params.panel_regions} | \
            awk -F '\\t' '$4 !~ /N/ && $5 !~ /N/' | \
            bcftools norm -m +any -O vcf -o {output.vcf}
        fi
        """

rule vcf2maf_annotate:
    input:
        vcf = rules.filter_and_restrict_sage.output.vcf
    output:
        vep_vcf = temp(os.path.join(SAGE_OUTDIR, "04-capturespace/{sample}.capspace.vep.vcf")),
        maf = os.path.join(SAGE_OUTDIR, "05-MAFs/{sample}.sage.maf")
    params:
        custom_enst = config["lcr-modules"]["cfDNA_SAGE_workflow"]["custom_enst"],
        vep_data = config["lcr-modules"]["cfDNA_SAGE_workflow"]["vep_data"],
        centre = config["lcr-modules"]["cfDNA_SAGE_workflow"]["centre"],
        normal_name = get_normal_name,
        ref_fasta = config["lcr-modules"]["_shared"]["ref_genome"],
        ref_ver = config["lcr-modules"]["_shared"]["ref_genome_ver"]
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
        --tumor-id {wildcards.sample} --normal-id {params.normal_name} \
        --vep-path $CONDA_PREFIX/bin/ \
        --vep-data {params.vep_data} --vep-forks {threads} \
        --custom-enst {params.custom_enst} --ref-fasta {params.ref_fasta} \
        --species homo_sapiens --ncbi-build {params.ref_ver} \
        --retain-info LPS,LPS_RC \
        --maf-center {params.centre} 2> {log} >> {log}
        """

rule recount_alleles:
    input:
        maf = rules.vcf2maf_annotate.output.maf,
        tbam = config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"]
    output:
        maf = temp(os.path.join(SAGE_OUTDIR, "08-augmentssm/{sample}.sage.augment.maf"))
    resources:
        mem_mb = 10000
    threads: 4
    conda:
        "envs/augment_ssm.yaml"
    params:
        script = os.path.join(UTILSDIR, "recount_alleles.py"),
        ref_genome_version = config["lcr-modules"]["_shared"]["ref_genome_ver"],
    shell: """
    python {params.script} \
    --bam {input.tbam} \
    --maf {input.maf} \
    --genome {params.ref_genome_version} \
    --output {output.maf} \
    --threads {threads}
    """

rule filter_repetitive_seq:
    """
    Remove mutations which are adjacent to a repetitive sequence.
    """
    input:
        maf = rules.recount_alleles.output.maf
    output:
        maf = temp(os.path.join(SAGE_OUTDIR, "10-filter_repeat/{sample}.sage.repeat_filt.maf"))
    params:
        max_repeat_len = 6,
        ref_fasta = config["lcr-modules"]["_shared"]["ref_genome"]
    conda:
        "envs/fastp.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        f"""python {os.path.join(UTILSDIR, "filter_rep_seq.py")} \
        --in_maf {{input.maf}} --out_maf {{output.maf}} --reference_genome {{params.ref_fasta}} --max_repeat_len {{params.max_repeat_len}}
        """

rule add_UMI_support:
    input:
        maf = rules.filter_repetitive_seq.output.maf,
        umi_bam = config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"]
    output:
        omaf = os.path.join(SAGE_OUTDIR, "11-umi/{sample}.umi_support.maf")
    params:
        script = os.path.join(UTILSDIR, "FetchVariantUMIs.py"),
        min_map_qual = config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_map_qual"]
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.add_UMI_support.log")
    shell:
        """python {params.script} --input_maf {input.maf} --input_bam {input.umi_bam} --output_maf {output.omaf} --min_map_quality {params.min_map_qual} &> {log}
        """

rule custom_filters:
    input:
        maf = rules.add_UMI_support.output.omaf
    output:
        maf = temp(os.path.join(SAGE_OUTDIR, "12-filtered/{sample}.sage.filtered.maf"))
    params:
        exac_freq = float(config["lcr-modules"]["cfDNA_SAGE_workflow"]["exac_max_freq"]),
        script = os.path.join(UTILSDIR, "custom_filters.py"),
        hotspot_vcf = config["lcr-modules"]["cfDNA_SAGE_workflow"]["hotspot_vcf"],
        blacklist_txt = config["lcr-modules"]["cfDNA_SAGE_workflow"]["blacklist_manifest"],
        min_germline_depth = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_germline_depth"]),
        min_alt_depth = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_alt_depth"]),
        min_tum_VAF = float(config["lcr-modules"]["cfDNA_SAGE_workflow"]["novel_vaf"]),
        min_t_depth = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_t_depth"]),
        min_UMI_3_count = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_UMI_3_count"]),
        low_alt_thresh = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["low_alt_thresh"]),
        low_alt_min_UMI_3_count = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["low_alt_min_UMI_3_count"]),
        background_rates = config["lcr-modules"]["cfDNA_SAGE_workflow"]["background_rates"],
        min_background_samples = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_background_samples"]),
        background_n_std = float(config["lcr-modules"]["cfDNA_SAGE_workflow"]["background_n_std"]),
        max_bg_rate = float(config["lcr-modules"]["cfDNA_SAGE_workflow"]["max_background_rate"])
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.custom_filters.log")
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        """
        python {params.script} --input_maf {input.maf} --output_maf {output.maf} \
        --min_tumour_vaf {params.min_tum_VAF} \
        --min_alt_depth_tum {params.min_alt_depth} \
        --min_germline_depth {params.min_germline_depth} \
        --min_t_depth {params.min_t_depth} \
        --gnomad_threshold {params.exac_freq} \
        --min_UMI_3_count {params.min_UMI_3_count} \
        --low_alt_thresh {params.low_alt_thresh} \
        --low_alt_min_UMI_3_count {params.low_alt_min_UMI_3_count} \
        $( [ -n "{params.blacklist_txt}" ] && echo "--blacklist {params.blacklist_txt}" || echo "" ) \
        $( [ -n "{params.hotspot_vcf}" ] && echo "--hotspots {params.hotspot_vcf}" || echo "" ) \
        $( [ -n "{params.background_rates}" ] && echo "--background_rates {params.background_rates}" || echo "" ) \
        --min_background_samples {params.min_background_samples} \
        --background_n_std {params.background_n_std} --max_background_rate {params.max_bg_rate} &> {log}
        """


rule augment_maf:
    input:
        index_maf = rules.custom_filters.output.maf,
        additional_mafs = older_sample_mafs,
        index_bam = config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"]
    output:
        maf = os.path.join(SAGE_OUTDIR, f"99-final/{{sample}}.v{pv}.tempest.maf")
    params:
        script = os.path.join(UTILSDIR, "augmentMAF.py"),
        ref_genome_version = config["lcr-modules"]["_shared"]["ref_genome_ver"],
        alt_support = config["lcr-modules"]["cfDNA_SAGE_workflow"]["min_alt_depth"],
        min_UMI_3_count = int(config["lcr-modules"]["cfDNA_SAGE_workflow"]["aug_min_UMI_3_count"]),
        phase_ID_col = config["lcr-modules"]["cfDNA_SAGE_workflow"]["phase_id_col"],
        phased_min_t_alt = config["lcr-modules"]["cfDNA_SAGE_workflow"]["phased_min_t_alt"]
    log:
        os.path.join(SAGE_OUTDIR, "logs/{sample}.augment_maf.log")
    conda:
        "envs/augment_ssm.yaml"
    resources:
        mem_mb = 20000
    threads: 3
    shell:
        """python {params.script} --sample_id {wildcards.sample} --index_maf {input.index_maf} --index_bam {input.index_bam} \
        --add_maf_files {input.additional_mafs} --genome_build {params.ref_genome_version} --threads {threads} \
        --min_alt_count {params.alt_support} \
        --phase_ID_col {params.phase_ID_col} --phased_min_t_alt_count {params.phased_min_t_alt} \
        --compute_umi_metrics --min_UMI_3_count {params.min_UMI_3_count} \
        --output {output.maf} &> {log}
        """

# Generate IGV screenshots of these variants
rule igv_screenshot_variants:
    input:
        maf = rules.augment_maf.output.maf,
        tbam = config["lcr-modules"]["cfDNA_SAGE_workflow"]["bam_path"],
        nbam = get_normal_bam
    output:
        html = os.path.join(SAGE_OUTDIR, "07-IGV/{sample}_report.html")
    params:
        refgenome = config["lcr-modules"]["_shared"]["ref_genome"],
        cytoband = config["lcr-modules"]["cfDNA_SAGE_workflow"]["cytoband"],
        genes = config["lcr-modules"]["cfDNA_SAGE_workflow"]["genetrack"],
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
        create_report --fasta {params.refgenome} --type mutation --tracks {input.tbam} {input.nbam} {params.genes} --flanking 1500 --output {output.html} --standalone --title {params.sample_name} --ideogram {params.cytoband} {input.maf} > {log}
        """

rule all_sage:
    input:
        expand(os.path.join(SAGE_OUTDIR, f"99-final/{{sample}}.v{pv}.tempest.maf"), sample=SAMPLESHEET.loc[SAMPLESHEET["tissue_status"]=="tumor"]["sample_id"]),
        expand(str(rules.igv_screenshot_variants.output.html), sample = SAMPLESHEET.loc[SAMPLESHEET["tissue_status"]=="tumor"]["sample_id"])