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
# `CFG` is a shortcut to `config["lcr-modules"]["deepsomatic"]`
CFG = op.setup_module(
    name = "deepsomatic",
    version = "1.0",
    subdirectories = ["inputs", "deepsomatic_temp", "deepsomatic_output", "filter", "gnomad", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _deepsomatic_input_bam,
    _deepsomatic_input_normal_bam,
    _deepsomatic_output_vcf,
    _cleanup_intermediate_dir,
    _deepsomatic_all


VERSION_MAP_DEEPSOMATIC = CFG["options"]["version_map"]
SEQTYPE_MAP_DEEPSOMATIC = CFG["options"]["seqtype_map"]

possible_genome_builds = ", ".join(list(VERSION_MAP_DEEPSOMATIC.keys()))
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        f"Samples table includes genome builds not yet compatible with this module. "
        f"This module is currently only compatible with {possible_genome_builds}. "
    )

possible_seq_types = ", ".join(list(SEQTYPE_MAP_DEEPSOMATIC.keys()))
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


# Input tumor bam
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _deepsomatic_input_bam:
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


# Input normal bam
# Symlinks the input files into the module results directory (under '00-inputs/')
rule _deepsomatic_input_normal_bam:
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


# Helper function to fetch normal_name for the unmatched normal bam file depending on sample chemistry
def get_normal(wildcards):
    CFG = config["lcr-modules"]["deepsomatic"]
    this_sample = op.filter_samples(
        CFG["runs"],
        sample_id = wildcards.tumour_id,
        seq_type = wildcards.seq_type, 
        genome_build = wildcards.genome_build
    )
    normal = this_sample["normal_name"].tolist()[0]
    return normal

# Fetch calling mode from config (either unmatched or tumor_only)
DEEPSOMATIC_CALLING_MODE = CFG["options"]["calling_mode"]

# Map calling mode to the appropriate DeepSomatic model
DEEPSOMATIC_MODEL_TYPES = {
    "unmatched": "ONT",
    "tumor_only": "ONT_TUMOR_ONLY",
}

# Map calling mode to the name used in output directories
DEEPSOMATIC_CALLING_MODE_DIRS = {
    "unmatched": "unmatched",
    "tumor_only": "tumour-only",
}

# Throw an error if calling_mode is not one of the accepted values
if DEEPSOMATIC_CALLING_MODE not in DEEPSOMATIC_MODEL_TYPES:
    raise ValueError(
        "Invalid DeepSomatic calling_mode: "
        f"{DEEPSOMATIC_CALLING_MODE}. "
        "Expected 'unmatched' or 'tumor_only'."
    )

DEEPSOMATIC_CALLING_MODE_DIR = DEEPSOMATIC_CALLING_MODE_DIRS[
    DEEPSOMATIC_CALLING_MODE
]


# Helper function to require the normal BAM only in unmatched mode
def _deepsomatic_get_normal_bam(wildcards):
    if DEEPSOMATIC_CALLING_MODE == "tumor_only":
        return []
    return str(rules._deepsomatic_input_normal_bam.output.bam)


# Helper function to supply normal-specific arguments only in unmatched mode
def _deepsomatic_get_normal_args(wildcards, input):
    if DEEPSOMATIC_CALLING_MODE == "tumor_only":
        return ""
    return (
        f'--reads_normal="{input.normal_bam}" '
        '--sample_name_normal="NORMAL"'
    )

# Supply filters requiring the normal if calling_mode is "tumor_only"
DEEPSOMATIC_FILTERS = CFG["options"]["filters"]
# Helper function to supply normal filters to filtering rule depending on calling_mode
if DEEPSOMATIC_CALLING_MODE == "tumor_only":
    DEEPSOMATIC_NORMAL_FILTER = ""
else:
    DEEPSOMATIC_NORMAL_FILTER = (
        f'&& FMT/NDP[0] >= {DEEPSOMATIC_FILTERS["normal_min_depth"]} '
        f'&& FMT/NAF[0:1] < {DEEPSOMATIC_FILTERS["normal_max_af"]}'
    )


# Call variants in tumour using DeepSomatic
rule _deepsomatic_call_variants:
    input:
        tumour_bam = str(rules._deepsomatic_input_bam.output.bam),
        normal_bam = _deepsomatic_get_normal_bam,
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        vcf = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.vcf.gz"
    log:
        stdout = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/deepsomatic.call_variants.stdout.log",
        stderr = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/deepsomatic.call_variants.stderr.log"
    params:
        model_type = DEEPSOMATIC_MODEL_TYPES[DEEPSOMATIC_CALLING_MODE],
        normal_args = _deepsomatic_get_normal_args,
        deepsomatic_args = CFG["options"]["deepsomatic_args"],
        intermediate_results_dir = CFG["dirs"]["deepsomatic_temp"] + "{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "--intermediate_results/",
        logging_dir = CFG["logs"]["deepsomatic_output"] + "logs/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR
    container:
        CFG["container_envs"]["deepsomatic"]
    threads:
        CFG["threads"]["deepsomatic"]
    resources:
        **CFG["resources"]["deepsomatic"]
    shell:
        op.as_one_line("""
        run_deepsomatic
            --model_type={params.model_type}
            --ref="{input.fasta}"
            --reads_tumor="{input.tumour_bam}"
            {params.normal_args}
            --output_vcf="{output.vcf}"
            --sample_name_tumor="TUMOR"
            --num_shards={threads}
            --logging_dir="{params.logging_dir}"
            --intermediate_results_dir="$(readlink -m "{params.intermediate_results_dir}")"
            {params.deepsomatic_args}
            > {log.stdout}
            2> {log.stderr}
        """)


# Create TBI index DeepSomatic VCF file
rule _deepsomatic_index:
    input:
        vcf = str(rules._deepsomatic_call_variants.output.vcf)
    output:
        tbi = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.vcf.gz.tbi"
    log:
        stdout = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/deepsomatic.index.stdout.log",
        stderr = CFG["dirs"]["deepsomatic_output"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/deepsomatic.index.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    shell:
        op.as_one_line("""
        tabix -p vcf {input.vcf} > {log.stdout} 2> {log.stderr}
        """)


# Filter out population variants and low quality variant calls
rule _deepsomatic_filter:
    input:
        vcf = str(rules._deepsomatic_call_variants.output.vcf), 
        tbi = str(rules._deepsomatic_index.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.filtered.vcf.gz",
        tbi = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.filtered.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/filter.stderr.log",
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/filter.stdout.log"
    params:
        min_depth = CFG["options"]["filters"]["min_depth"],
        snv_min_vaf = CFG["options"]["filters"]["snv_min_vaf"],
        snv_min_alt_depth = CFG["options"]["filters"]["snv_min_alt_depth"],
        indel_min_vaf = CFG["options"]["filters"]["indel_min_vaf"],
        indel_min_alt_depth = CFG["options"]["filters"]["indel_min_alt_depth"],
        normal_filter = DEEPSOMATIC_NORMAL_FILTER
    shell:
        op.as_one_line("""
        bcftools isec -C -w1 {input.vcf} {input.pon} 2> {log.stderr} |
        bcftools view
            -i 'FILTER="PASS" &&
                FMT/DP[0] >= {params.min_depth} &&
                (
                    (TYPE="snp" && FMT/VAF[0:0] >= {params.snv_min_vaf} && FMT/AD[0:1] >= {params.snv_min_alt_depth}) ||
                    (TYPE="indel" && FMT/VAF[0:0] >= {params.indel_min_vaf} && FMT/AD[0:1] >= {params.indel_min_alt_depth})
                )
                {params.normal_filter}'
            -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Annotates VCF file with gnomAD frequency data
rule _deepsomatic_gnomad_annotation:
    input:
        vcf = str(rules._deepsomatic_filter.output.vcf),
        tbi = str(rules._deepsomatic_filter.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.gnomad.vcf.gz",
        tbi = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/output.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/gnomad.stderr.log",
        stdout = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + "/gnomad.stdout.log"
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
        awk 'BEGIN {{FS=OFS="\t"}} /^#/ {{print; next}} {{ if ($8 == ".") $8="AF=0"; else if ($8 !~ /(^|;)AF=/) $8=$8";AF=0"; print }}' |
        bcftools view -i 'INFO/AF[0] < 0.0001' -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _deepsomatic_output_vcf:
    input:
        vcf = str(rules._deepsomatic_gnomad_annotation.output.vcf),
        tbi = str(rules._deepsomatic_gnomad_annotation.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + ".output.final.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR + ".output.final.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Remove files from the intermediate_results_dir used during DeepSomatic
rule _cleanup_intermediate_dir:
    input:
        vcf = str(rules._deepsomatic_output_vcf.output.vcf),
        tbi = str(rules._deepsomatic_output_vcf.output.tbi)
    output:
        cleanup_dummy = touch(
            CFG["dirs"]["deepsomatic_output"]
            + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}"
            + "--{chemistry}--" + DEEPSOMATIC_CALLING_MODE_DIR
            + "/cleanup_complete.txt"
        )
    params:
        cleanup_toggle = CFG["options"]["cleanup_toggle"],
        intermediate_results_parent = CFG["dirs"]["deepsomatic_temp"],
        intermediate_results_dir = (
            CFG["dirs"]["deepsomatic_temp"]
            + "{tumour_id}--{normal_name}--{chemistry}--"
            + DEEPSOMATIC_CALLING_MODE_DIR
            + "--intermediate_results/"
        )
    shell:
        op.as_one_line("""
        if [[ "{params.cleanup_toggle}" == "True" ]]; then
            target=$(readlink -f "{params.intermediate_results_dir}" 2>/dev/null || true) ;
            parent=$(readlink -f "{params.intermediate_results_parent}" 2>/dev/null || true) ;
            if [[ -n "$target" &&
                  -n "$parent" &&
                  "$target" != "/" &&
                  "$target" == "$parent"/* ]]; then
                rm -rf -- "$target" ;
            elif [[ -n "$target" ]]; then
                echo "Refusing to remove unexpected path: $target" >&2 ;
                exit 1 ;
            fi ;
            if [[ -L "{params.intermediate_results_dir}" ]]; then
                rm -f -- "{params.intermediate_results_dir}" ;
            fi ;
        else
            echo "Skipping cleanup" ;
        fi
        """)


# Generates the target sentinels for each run, which generate the symlinks
rule _deepsomatic_all:
    input:
        expand(
            [
                str(rules._deepsomatic_output_vcf.output.vcf),
                str(rules._deepsomatic_output_vcf.output.tbi),
                str(rules._cleanup_intermediate_dir.output.cleanup_dummy)
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            chemistry=CFG["runs"]["tumour_chemistry"],
            normal_name=CFG["runs"]["normal_name"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
