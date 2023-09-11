#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  n/a
# Module Author:    Krysta M Coyle PhD
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op

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
# `CFG` is a shortcut to `config["lcr-modules"]["sv_repair"]`
CFG = op.setup_module(
    name = "sv_repair",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs","filter_vcf", "vcf_to_df","bed","rmsk","nonb_dna", "repair", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _sv_repair_input_vcf,
    _sv_repair_filter_vcf,
    #_sv_repair_output_unk


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')

rule _sv_repair_input_vcf:
    input:
        vcf = CFG["inputs"]["sample_vcf"]
    output:
        vcf = CFG["dirs"]["inputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.vcf"
    run:
        op.absolute_symlink(input.vcf, output.vcf)


rule _sv_repair_filter_vcf:
    input:
        vcf = str(rules._sv_repair_input_vcf.output.vcf)
    output:
        file = CFG["dirs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter.complete",
        vcf = CFG["dirs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise.recode.vcf"
    log:
        stdout = CFG["logs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/filter_vcf.stdout.log"
    params:
        prefix = CFG["dirs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise"
    conda:
        CFG["conda_envs"]["vcftools"]
    shell:
        op.as_one_line("""
            vcftools --vcf {input.vcf} --remove-filtered-all --remove-INFO "IMPRECISE" --recode --recode-INFO-all --out {params.prefix} && touch {output.file}
        """)

rule _sv_repair_vcf_to_df:
    input:
        vcf = CFG["dirs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise.recode.vcf"
    output:
        tsv = CFG["dirs"]["vcf_to_df"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise.tsv"
    log:
        stdout = CFG["logs"]["vcf_to_df"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vcf_to_df.stdout.log"
    script:
        "src/python/vcf.fields.py"

rule _sv_repair_vcf_to_bed:
    input:
        vcf = CFG["dirs"]["filter_vcf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise.recode.vcf"
    output:
        bed = CFG["dirs"]["bed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/somaticSV.pass.precise.bed"
    log:
        stdout = CFG["logs"]["bed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/vcf_to_bed.stdout.log"
    shell:
        op.as_one_line("""
            awk '{{if ($0 ~ /^[[:space:]]*#/){{NR--}}else{{print $1"\t"$2-1"\t"$2}}}}' {input.vcf} > {output.bed}
        """)

rule _sv_repair_bed_rmsk:
    input:
        bed = str(rules._sv_repair_vcf_to_bed.output.bed),
        tsv = str(rules._sv_repair_vcf_to_df.output.tsv),
        rmsk = reference_files("genomes/{genome_build}/repeatmasker/repeatmasker.{genome_build}.bed")
    output:
        bed = CFG["dirs"]["rmsk"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/intersect.rmsk.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    log:
        stdout = CFG["logs"]["rmsk"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/intersect_rmsk.stdout.log"
    shell:
        op.as_one_line("""
            bedtools intersect -a {input.bed} -b {input.rmsk} -loj > {output.bed}
        """)

rule _sv_repair_nonb_dna:
    input:
        ref_bed = reference_files("genomes/{genome_build}/nonb_dna/{genome_build}.bed"),
        bed = str(rules._sv_repair_vcf_to_bed.output.bed),
        tsv = str(rules._sv_repair_vcf_to_df.output.tsv)
    output:
        bed = CFG["dirs"]["nonb_dna"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/intersect.nonb.bed"
    conda:
        CFG["conda_envs"]["bedtools"]
    resources:
        mem_mb = 10000
    log:
        stdout = CFG["logs"]["nonb_dna"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/intersect_nonb.stdout.log"
    shell:
        op.as_one_line("""
            bedtools intersect -a {input.bed} -b {input.ref_bed} -loj > {output.bed}
        """)

rule _sv_repair_mechanisms:
    input:
        vcf_df = str(rules._sv_repair_vcf_to_df.output.tsv),
        rmsk = str(rules._sv_repair_bed_rmsk.output.bed),
        non_b = str(rules._sv_repair_nonb_dna.output.bed)
    params:
        name = "{tumour_id}--{normal_id}--{pair_status}"
    output:
        tsv = CFG["dirs"]["repair"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/repair.mech.tsv",
        paired = CFG["dirs"]["repair"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/repair.table.tsv"
    resources:
        mem_mb = 10000
    log:
        stdout = CFG["logs"]["repair"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/repair.mech.log"
    script:
        "src/R/repair_mech.R"

rule _sv_output_files:
    input:
        tsv = str(rules._sv_repair_mechanisms.output.tsv),
        paired = str(rules._sv_repair_mechanisms.output.paired)
    output:
        tsv = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/repair.mech.tsv",
        paired = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/repair.table.tsv"
    run:
        op.absolute_symlink(input.tsv, output.tsv),
        op.absolute_symlink(input.paired, output.paired)

# Generates the target sentinels for each run, which generate the symlinks
rule _sv_repair_all:
    input:
        expand(
            [
                str(rules._sv_output_files.output.paired),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"]),
        expand(
            [
                str(rules._sv_output_files.output.tsv),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
