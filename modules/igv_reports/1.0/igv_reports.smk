#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


import re
import oncopipe as op

min_oncopipe_version = "1.0.11"
import pkg_resources
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

current_version = pkg_resources.get_distribution("oncopipe").version
if version.parse(current_version) < version.parse(min_oncopipe_version):
    logger.warning(
        '\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}'
        "\n" f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m'
    )
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

CFG = op.setup_module(
    name = "igv_reports",
    version = "1.0",
    subdirectories = ["inputs", "filtered_maf", "igv_reports", "outputs"],
)

# Pre-capture values needed in lambdas — CFG is deleted by op.cleanup_module
# and would be out of scope when Snakemake evaluates lambdas at DAG build time.
_igv_tools = list(CFG["inputs"]["mafs"].keys())
_igv_mafs = dict(CFG["inputs"]["mafs"])
_igv_info_columns = " ".join(CFG["options"]["info_columns"])
_igv_genome_map = dict(CFG["options"]["genome_map"])

localrules:
    _igv_reports_input_maf,
    _igv_reports_input_bam,
    _igv_reports_output_html,
    _igv_reports_all


##### RULES #####


rule _igv_reports_input_maf:
    input:
        maf = lambda w: _igv_mafs[w.tool].format(
            seq_type=w.seq_type, genome_build=w.genome_build,
            tumour_id=w.tumour_id, normal_id=w.normal_id,
            pair_status=w.pair_status
        )
    output:
        maf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.maf"
    wildcard_constraints:
        tool = "|".join(re.escape(t) for t in _igv_tools)
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _igv_reports_input_bam:
    input:
        tumour_bam = CFG["inputs"]["tumour_bam"],
        normal_bam = CFG["inputs"]["normal_bam"],
    output:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.tumour.bam",
        tumour_bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.tumour.bam.bai",
        tumour_crai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.tumour.bam.crai",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.normal.bam",
        normal_bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.normal.bam.bai",
        normal_crai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}.normal.bam.crai",
    run:
        op.absolute_symlink(input.tumour_bam, output.tumour_bam)
        op.absolute_symlink(input.tumour_bam + ".bai", output.tumour_bai)
        op.absolute_symlink(input.tumour_bam + ".bai", output.tumour_crai)
        op.absolute_symlink(input.normal_bam, output.normal_bam)
        op.absolute_symlink(input.normal_bam + ".bai", output.normal_bai)
        op.absolute_symlink(input.normal_bam + ".bai", output.normal_crai)


rule _igv_reports_filter_maf:
    input:
        maf = str(rules._igv_reports_input_maf.output.maf)
    output:
        maf = CFG["dirs"]["filtered_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.drivers.maf"
    log:
        CFG["logs"]["filtered_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}/filter_maf.log"
    conda:
        CFG["conda_envs"]["pandas"]
    threads:
        CFG["threads"]["filter_maf"]
    resources:
        **CFG["resources"]["filter_maf"]
    wildcard_constraints:
        tool = "|".join(re.escape(t) for t in _igv_tools)
    params:
        driver_genes = CFG["options"]["driver_genes"],
        driver_gene_col = CFG["options"]["driver_gene_col"],
        driver_col = CFG["options"]["driver_col"],
        driver_col_value = CFG["options"]["driver_col_value"],
        maf_gene_col = CFG["options"]["maf_gene_col"],
        maf_vc_col = CFG["options"]["maf_vc_col"],
        coding_variant_classes = CFG["options"]["coding_variant_classes"]
    script:
        "src/filter_maf_for_igv.py"


rule _igv_reports_run:
    input:
        maf = str(rules._igv_reports_filter_maf.output.maf),
        tumour_bam = str(rules._igv_reports_input_bam.output.tumour_bam),
        tumour_bai = str(rules._igv_reports_input_bam.output.tumour_bai),
        tumour_crai = str(rules._igv_reports_input_bam.output.tumour_crai),
        normal_bam = str(rules._igv_reports_input_bam.output.normal_bam),
        normal_bai = str(rules._igv_reports_input_bam.output.normal_bai),
        normal_crai = str(rules._igv_reports_input_bam.output.normal_crai)
    output:
        html = CFG["dirs"]["igv_reports"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.html"
    log:
        CFG["logs"]["igv_reports"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}/igv_reports.log"
    conda:
        CFG["conda_envs"]["igv_reports"]
    threads:
        CFG["threads"]["igv_reports"]
    resources:
        **CFG["resources"]["igv_reports"]
    wildcard_constraints:
        tool = "|".join(re.escape(t) for t in _igv_tools)
    params:
        flanking = CFG["options"]["flanking"],
        info_columns = _igv_info_columns,
        title = lambda w: f"{w.tumour_id} vs {w.normal_id} — {w.tool} drivers",
        igv_genome = lambda w: _igv_genome_map.get(w.genome_build, w.genome_build)
    shell:
        op.as_one_line("""
        create_report {input.maf}
            --genome {params.igv_genome}
            --tracks {input.tumour_bam} {input.normal_bam}
            --flanking {params.flanking}
            --info-columns {params.info_columns}
            --title "{params.title}"
            --output {output.html}
            > {log} 2>&1
        """)


rule _igv_reports_output_html:
    input:
        html = str(rules._igv_reports_run.output.html)
    output:
        html = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{tool}.html"
    wildcard_constraints:
        tool = "|".join(re.escape(t) for t in _igv_tools)
    run:
        op.relative_symlink(input.html, output.html, in_module=True)


rule _igv_reports_all:
    input:
        expand(
            expand(
                str(rules._igv_reports_output_html.output.html),
                tool=_igv_tools,
                allow_missing=True
            ),
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        )


##### CLEANUP #####


op.cleanup_module(CFG)
