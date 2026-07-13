#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giuliano Banco
# Module Author:    Giuliano Banco
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
import os

# Check that the oncopipe dependency is up-to-date. Add all the following lines to any module that uses new features in oncopipe
min_oncopipe_version="1.0.11"
from importlib.metadata import version as pkg_version
try:
    from packaging import version
except ModuleNotFoundError:
    sys.exit("The packaging module dependency is missing. Please install it ('pip install packaging') and ensure you are using the most up-to-date oncopipe version")

# To avoid this we need to add the "packaging" module as a dependency for LCR-modules or oncopipe

current_version = pkg_version("oncopipe")
if version.parse(current_version) < version.parse(min_oncopipe_version):
    print('\x1b[0;31;40m' + f'ERROR: oncopipe version installed: {current_version}' + '\x1b[0m')
    print('\x1b[0;31;40m' + f"ERROR: This module requires oncopipe version >= {min_oncopipe_version}. Please update oncopipe in your environment" + '\x1b[0m')
    sys.exit("Instructions for updating to the current version of oncopipe are available at https://lcr-modules.readthedocs.io/en/latest/ (use option 2)")

# End of dependency checking section

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["clairs"]`
CFG = op.setup_module(
    name = "clairs",
    version = "1.1",
    subdirectories = ["inputs", "clairs", "filter", "gnomad", "outputs"]
)

config["pipeline_name"] = "clairs.yaml"

# Define rules to be run locally when using a compute cluster
localrules:
    _clairs_input_bam,
    _clairs_input_normal_bam,
    _clairs_output_vcf,
    _clairs_clean,
    _clairs_all


VERSION_MAP_CLAIRS = CFG["options"]["version_map"]
SEQTYPE_MAP_CLAIRS = CFG["options"]["seqtype_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP_CLAIRS.keys()))
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        f"Samples table includes genome builds not yet compatible with this module. "
        f"This module is currently only compatible with {possible_genome_builds}. "
    )

possible_seq_types = ", ".join(list(SEQTYPE_MAP_CLAIRS.keys()))
for seq_type in CFG["runs"]["tumour_seq_type"]:
    assert seq_type in possible_seq_types, (
        f"Samples table includes seq types not yet compatible with this module. "
        f"This module is currently only compatible with {possible_seq_types}. "
    )


chemistry = CFG["runs"]["tumour_chemistry"]

normal_map = CFG["options"]["normal_name"]

CFG["runs"]["normal_name"] = CFG["runs"]["tumour_chemistry"].map(normal_map)

if CFG["runs"]["normal_name"].isna().any():
    missing = CFG["runs"].loc[CFG["runs"]["normal_name"].isna(), "tumour_chemistry"].unique()
    print(f"Warning: some chemistry values have no mapping in normal_name. Available chemistry values: R9, R10")

##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _clairs_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


rule _clairs_input_normal_bam:
    input:
        bam = CFG["inputs"]["normal_bam"],
        bai = CFG["inputs"]["normal_bai"],
        normal_name = CFG["runs"]["normal_name"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_name}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_name}.bam.bai",
        crai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_name}.bam.crai"
    run:
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)
        op.absolute_symlink(input.bai, output.crai)


def get_platform(wildcards):
    CFG = config["lcr-modules"]["clairs"]
    this_sample = op.filter_samples(
        CFG["runs"],
        sample_id = wildcards.tumour_id,
        seq_type = wildcards.seq_type, 
        genome_build = wildcards.genome_build
    )
    platform = this_sample["tumour_platform"].tolist()[0]
    return platform


def get_normal(wildcards):
    CFG = config["lcr-modules"]["clairs"]
    this_sample = op.filter_samples(
        CFG["runs"],
        sample_id = wildcards.tumour_id,
        seq_type = wildcards.seq_type, 
        genome_build = wildcards.genome_build
    )
    normal = this_sample["normal_name"].tolist()[0]
    return normal


# Calls variants using ClairS
rule _clairs_call_variants:
    input:
        tumour_bam = str(rules._clairs_input_bam.output.bam),
        normal_bam = str(rules._clairs_input_normal_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
    output:
        vcf = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.vcf.gz",
        tbi = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.stdout.log",
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.stderr.log"
    params:
        clairs_args = CFG["options"]["clairs_args"],
        platform = get_platform,
        output_dir = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/",
        call_indels = lambda wc: (
            "--enable_indel_calling"
            if config["lcr-modules"]["clairs"]["options"].get("call_indels", False) else ""
            )
    container:
        CFG["container_envs"]["clairs"]
    threads:
        CFG["threads"]["clairs"]
    resources:
        **CFG["resources"]["clairs"]
    shell:
        op.as_one_line("""
        /opt/bin/run_clairs
            --tumor_bam_fn {input.tumour_bam}
            --normal_bam_fn {input.normal_bam}
            --ref_fn {input.fasta}
            --threads {threads}
            --platform {params.platform}
            --output_dir {params.output_dir}
            --conda_prefix /opt/conda/envs/clairs
            -s TUMOR
            {params.call_indels}
            {params.clairs_args}
            > {log.stdout} 2> {log.stderr}
        """)


# Combines ClairS VCF output files so indels and SNVs are in the same file
rule _clairs_combine_vcfs:
    input:
        vcf = str(rules._clairs_call_variants.output.vcf),
        tbi = str(rules._clairs_call_variants.output.tbi)
    output:
        vcf = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combined.vcf.gz",
        tbi = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combined.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combine_vcfs.stdout.log",
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combine_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    params:
        output_dir = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/",
        indels_enabled = lambda wc: (
            "enabled"
            if config["lcr-modules"]["clairs"]["options"].get("call_indels", False) else ""
            )
    shell:
        op.as_one_line("""
        if [ -z "{params.indels_enabled}" ]; then
            bcftools view -h {input.vcf} > {params.output_dir}/indel.vcf &&
            bgzip -f {params.output_dir}/indel.vcf &&
            tabix -p vcf {params.output_dir}/indel.vcf.gz ;
        fi ;
        bcftools sort
            {params.output_dir}/output.vcf.gz
            -Oz
            -o {params.output_dir}/output.sorted.vcf.gz
            >> {log.stdout} 2>> {log.stderr}
        &&
        tabix -p vcf
            {params.output_dir}/output.sorted.vcf.gz
            >> {log.stdout} 2>> {log.stderr}
        &&
        bcftools sort
            {params.output_dir}/indel.vcf.gz
            -Oz
            -o {params.output_dir}/indel.sorted.vcf.gz
            >> {log.stdout} 2>> {log.stderr}
        &&
        tabix -p vcf
            {params.output_dir}/indel.sorted.vcf.gz
            >> {log.stdout} 2>> {log.stderr}
        &&
        bcftools concat
            --allow-overlaps
            {params.output_dir}/output.sorted.vcf.gz
            {params.output_dir}/indel.sorted.vcf.gz
            -Ou
            2>> {log.stderr}
        |
        bcftools sort
            -Oz
            -o {output.vcf}
            >> {log.stdout} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf}
            >> {log.stdout} 2>> {log.stderr}
        """)


# Cleans up additional files created by ClairS
rule _clairs_clean:
    input:
        str(rules._clairs_call_variants.output.vcf),
        str(rules._clairs_call_variants.output.tbi),
        str(rules._clairs_combine_vcfs.output.vcf),
        str(rules._clairs_combine_vcfs.output.tbi)
    output:
        cleanup = touch(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/cleanup_complete.txt")
    params:
        cleanup_toggle = CFG["options"]["cleanup_toggle"]
    shell:
        op.as_one_line("""
        d=$(dirname {output.cleanup}) &&
        if [ "{params.cleanup_toggle}" = "True" ] || [ "{params.cleanup_toggle}" = "true" ]; then
            rm -rf "$d/tmp/" &&
            rm -rf "$d/tmp_TUMOR/" &&
            rm -rf "$d/logs/" &&
            rm -f "$d"/run_clairs.log* &&
            rm -f "$d"/output.sorted.vcf.gz* &&
            rm -f "$d"/indel.sorted.vcf.gz* ;
        else
            echo "cleanup_toggle is false; Skipping cleanup" >&2 ;
        fi
        """)


# Filters out variants
rule _clairs_filter:
    input:
        vcf = str(rules._clairs_combine_vcfs.output.vcf), 
        tbi = str(rules._clairs_combine_vcfs.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.final.vcf.gz",
        tbi = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.final.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/filter_pon.stderr.log",
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/filter_pon.stdout.log"
    params:
        min_depth = CFG["options"]["filters"]["tumour_min_depth"],
        snv_min_af = CFG["options"]["filters"]["snv_min_af"],
        snv_min_alt_depth = CFG["options"]["filters"]["snv_min_alt_depth"],
        indel_min_af = CFG["options"]["filters"]["indel_min_af"],
        indel_min_alt_depth = CFG["options"]["filters"]["indel_min_alt_depth"],
        normal_min_depth = CFG["options"]["filters"]["normal_min_depth"],
        normal_max_af = CFG["options"]["filters"]["normal_max_af"]
    shell:
        op.as_one_line("""
        bcftools isec -C -w1 {input.vcf} {input.pon} 2> {log.stderr} |
        bcftools view
            -i 'FILTER="PASS" &&
                FMT/DP[0] >= {params.min_depth} &&
                (
                    (TYPE="snp" && FMT/AF[0] >= {params.snv_min_af} && FMT/AD[0:1] >= {params.snv_min_alt_depth}) ||
                    (TYPE="indel" && FMT/AF[0] >= {params.indel_min_af} && FMT/AD[0:1] >= {params.indel_min_alt_depth})
                ) &&
                FMT/NDP[0] >= {params.normal_min_depth} &&
                FMT/NAF[0] < {params.normal_max_af}'
            -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Annotates VCF file with gnomAD frequency data and filters out poor calls
rule _clairs_gnomad_annotation:
    input:
        vcf = str(rules._clairs_filter.output.vcf),
        tbi = str(rules._clairs_filter.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.gnomad.vcf.gz",
        tbi = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs_gnomad_annotation.stderr.log",
        stdout = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs_gnomad_annotation.stdout.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads}
        -a {input.gnomad} -c INFO/AF {input.vcf} 2> {log.stderr} |
        awk 'BEGIN {{FS=OFS="\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' |
        bcftools view -i 'INFO/AF < 0.0001' -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _clairs_output_vcf:
    input:
        vcf = str(rules._clairs_gnomad_annotation.output.vcf),
        tbi = str(rules._clairs_gnomad_annotation.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched.clairs.combined.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched.clairs.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _clairs_all:
    input:
        expand(
            [
                str(rules._clairs_output_vcf.output.vcf),
                str(rules._clairs_output_vcf.output.tbi),
                str(rules._clairs_clean.output.cleanup)
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            platform=CFG["runs"]["tumour_platform"],
            chemistry=CFG["runs"]["tumour_chemistry"],
            normal_name=CFG["runs"]["normal_name"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
