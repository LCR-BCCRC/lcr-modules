import pandas as pd
import os
import datetime
import sys
MODULE_PATH = os.path.join(config["lcr-modules"]["_shared"]["lcr-modules"], "modules/cfdna_pipeline/1.0/")
sys.path.append(MODULE_PATH) # add local module to path


TODAY = datetime.datetime.now().strftime("%m/%d/%Y")
BAM_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "bam_pipeline")
UTILSDIR = os.path.join(MODULE_PATH, "utils")
SAGE_OUTDIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "sage_pipeline")

COMPILE_REPORT_SCRIPT = os.path.join(MODULE_PATH, "patient_reports/compile_report.py")
REPORT_TEMP = config["lcr-modules"]["cfDNA_patient_reports"]["report_template"]
REPORTS_DIR = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "reports")
all_samples = config["lcr-modules"]["_shared"]["samples"].copy()
REP_SAMPLESHEET = all_samples.loc[all_samples["matched_normal"] != "unmatched"].copy()
# get list of patient_ids for patients who have a tumor sample
PATS_TO_REPORT = REP_SAMPLESHEET[REP_SAMPLESHEET["tissue_status"]== "tumor"]["patient_id"].unique().tolist()

localrules:
    record_sample_completion,
    _vcf2maf_crossmap

# input functions
def find_sage_outputs(wildcards):
    # make list of all sample names belonging to patient 
    patient_samples = REP_SAMPLESHEET[(REP_SAMPLESHEET["patient_id"] == wildcards.patient) & (REP_SAMPLESHEET['timepoint'] != 'normal' )]["sample_id"].tolist()
    return expand(os.path.join(SAGE_OUTDIR,"99-final/{sample}.processed.maf"), sample=patient_samples)

def find_completion_time(wildcards):
    patient_samples = REP_SAMPLESHEET[(REP_SAMPLESHEET["patient_id"] == wildcards.patient) & (REP_SAMPLESHEET['timepoint'] != 'normal' )]["sample_id"].tolist()
    return expand(os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "completion", "{sample}.completion.txt"), sample=patient_samples )

def find_hsmetrics(wildcards):
    patient_samples = REP_SAMPLESHEET[REP_SAMPLESHEET["patient_id"] == wildcards.patient]["sample_id"].unique().tolist()
    return expand(os.path.join(BAM_OUTDIR, "Q2-hs_metrics/{sample}.hs_metrics.txt"), sample=patient_samples)

def find_targ_cov(wildcards):
    # make list of all sample names belonging to patient
    patient_samples = REP_SAMPLESHEET[REP_SAMPLESHEET["patient_id"] == wildcards.patient]["sample_id"].unique().tolist()
    return expand(os.path.join(BAM_OUTDIR , "Q2-hs_metrics" , "{sample}.target_coverage.txt"), sample=patient_samples)

def find_igv_report(wildcards):
    patient_samples = REP_SAMPLESHEET[REP_SAMPLESHEET["patient_id"] == wildcards.patient]["sample_id"].unique().tolist()
    return expand(os.path.join(SAGE_OUTDIR, "07-IGV/{sample}_report.html"), sample=patient_samples)

def find_insert_length(wildcards):
    patient_samples = REP_SAMPLESHEET[REP_SAMPLESHEET["patient_id"] == wildcards.patient]["sample_id"].unique().tolist()
    return expand(os.path.join(BAM_OUTDIR, "Q4-insert_size", "{sample}.insert_size_metrics.txt"), sample=patient_samples)

def lymphgen_outputs(wildcards):
    """input function for lymphgen module, calls the outpts
    for _lymphgen_output_txt rule"""
    patient_samples = REP_SAMPLESHEET[(REP_SAMPLESHEET["patient_id"] == wildcards.patient) & (REP_SAMPLESHEET["tissue_status"] != "normal")]["sample_id"].unique().tolist()

    # if lymphgen rule is in workflow, then return the expand of hte files
    if config["lcr-modules"]["cfDNA_patient_reports"]["include_lymphgen"]:
        return expand(expand(str(rules._lymphgen_output_txt.output),zip,
                    seq_type="capture",
                    genome_build="hg38",
                    normal_id="None",
                    pair_status="no_normal",
                    cnvs_wc="no_cnvs",
                    sv_wc="no_sv",
                    A53_wc="no_A53", allow_missing=True),
                    tumour_id= patient_samples,)
    else:
        # make empty dumby file to return to rule
        empty_file = os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "empty_lymphen.txt")
        with open(empty_file, "w") as f:
            f.write("empty")
        return empty_file

        

rule record_sample_completion:
    output:
        os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "completion", "{sample}.completion.txt")
    shell:
        f"""
        # if file exists, touch it, so it is never overwritten
        if [ ! -f {{output}} ]; then
            echo "{TODAY}" > {{output}}
        else
            touch {{output}}
        fi
        """

rule compile_report:
    input:
        hs_metrics = find_hsmetrics,
        targ_cov = find_targ_cov,
        sage_calls = find_sage_outputs,
        sample_comp = find_completion_time,
        lg_status = lymphgen_outputs
    output:
        compiled_nb = os.path.join(REPORTS_DIR, "{patient}", "{patient}_report.ipynb")
    params:
        samplesheet_path = config["lcr-modules"]["_shared"]["samplesheet"]
    resources:
        mem_mb = 5000
    threads: 1
    conda:
        "envs/quarto.yaml"
    log:
        os.path.join(REPORTS_DIR, "logs" , "cmp_rep_{patient}.log")
    shell:
        f"""python {COMPILE_REPORT_SCRIPT} --in_notebook {REPORT_TEMP} --out_notebook {{output.compiled_nb}} --completion_files {{input.sample_comp}} \
        --maf_files {{input.sage_calls}} --samplesheet_path {{params.samplesheet_path}} --repo_path {MODULE_PATH} --lymphgen_output {{input.lg_status}} \
        --patient_id {{wildcards.patient}} --hs_metrics {{input.hs_metrics}} --targ_cov {{input.targ_cov}} &> {{log}}"""

rule convert_report_to_html:
    input:
        rules.compile_report.output.compiled_nb
    output:
        html_report = os.path.join(REPORTS_DIR, "{patient}/{patient}_report.html")
    conda:
        "envs/quarto.yaml"
    resources:
        mem_mb = 5000
    threads: 1
    shell:
        f"""
        cd {REPORTS_DIR}/{{wildcards.patient}}
        quarto render {{wildcards.patient}}_report.ipynb -o {{wildcards.patient}}_report.html --to html --quiet
        """

# need to include cause lymphgen runs on grch37
rule _vcf2maf_crossmap:
    input:
        maf = os.path.join(SAGE_OUTDIR, "99-final/{tumour_id}.processed.maf"),
        convert_coord = config['lcr-modules']["lymphgen"]["convert_coord"],
        chains = config["lcr-modules"]["lymphgen"]["hg38_chainfile"]
    output:
        maf = os.path.join(SAGE_OUTDIR, "99-final/{tumour_id}.grch37.processed.maf"),
        bed = temp(os.path.join(SAGE_OUTDIR, "99-final/{tumour_id}.grch37.processed.unmapped.bed"))
    log:
        os.path.join(config["lcr-modules"]["_shared"]["root_output_dir"], "crossmap/logs" , "crossmap_{tumour_id}.log")
    conda:
        "envs/crossmap.yaml"
    threads: 1
    resources:
        mem_mb= 8000
    shell:
        """
        {input.convert_coord} \
        {input.maf} \
        {input.chains} \
        {output.maf} \
        crossmap &> {log}
        """

# add rule all to call outputs
rule make_all_reports:
    input:
        expand(str(rules.convert_report_to_html.output.html_report), patient =PATS_TO_REPORT)
