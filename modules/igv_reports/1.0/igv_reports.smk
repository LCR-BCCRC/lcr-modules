#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Ryan Morin
# Module Author:    Ryan Morin
# Contributors:     N/A


##### SETUP #####


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

localrules:
    _igv_reports_input_maf,
    _igv_reports_input_bam,
    _igv_reports_output_html,
    _igv_reports_all


##### RULES #####


rule _igv_reports_input_maf:
    input:
        maf = CFG["inputs"]["maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


rule _igv_reports_input_bam:
    input:
        tumour_bam = CFG["inputs"]["tumour_bam"],
        tumour_bai = CFG["inputs"]["tumour_bai"],
        normal_bam = CFG["inputs"]["normal_bam"],
        normal_bai = CFG["inputs"]["normal_bai"],
    output:
        tumour_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam",
        tumour_bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}.bam.bai",
        normal_bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam",
        normal_bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{normal_id}.bam.bai",
    run:
        op.absolute_symlink(input.tumour_bam, output.tumour_bam)
        op.absolute_symlink(input.tumour_bai, output.tumour_bai)
        op.absolute_symlink(input.normal_bam, output.normal_bam)
        op.absolute_symlink(input.normal_bai, output.normal_bai)


rule _igv_reports_filter_maf:
    input:
        maf = str(rules._igv_reports_input_maf.output.maf)
    output:
        maf = CFG["dirs"]["filtered_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.drivers.maf"
    log:
        CFG["logs"]["filtered_maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_maf.log"
    conda:
        CFG["conda_envs"]["pandas"]
    threads:
        CFG["threads"]["filter_maf"]
    resources:
        **CFG["resources"]["filter_maf"]
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
        normal_bam = str(rules._igv_reports_input_bam.output.normal_bam),
        normal_bai = str(rules._igv_reports_input_bam.output.normal_bai),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        html = CFG["dirs"]["igv_reports"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.html"
    log:
        CFG["logs"]["igv_reports"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/igv_reports.log"
    conda:
        CFG["conda_envs"]["igv_reports"]
    threads:
        CFG["threads"]["igv_reports"]
    resources:
        **CFG["resources"]["igv_reports"]
    params:
        flanking = CFG["options"]["flanking"],
        info_columns = lambda w: " ".join(CFG["options"]["info_columns"]),
        title = lambda w: f"{w.tumour_id} vs {w.normal_id}"
    shell:
        op.as_one_line("""
        create_report {input.maf}
            --fasta {input.fasta}
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
        html = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.html"
    run:
        op.relative_symlink(input.html, output.html, in_module=True)


rule _igv_reports_all:
    input:
        expand(
            [str(rules._igv_reports_output_html.output.html)],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]
        )


##### CLEANUP #####


op.cleanup_module(CFG)
