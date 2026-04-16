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
# `CFG` is a shortcut to `config["lcr-modules"]["clairs"]`
CFG = op.setup_module(
    name = "clairs",
    version = "1.0",
    subdirectories = ["inputs", "clairs", "gnomad", "filter", "outputs"],
)

config["pipeline_name"] = "clairs.yaml"

# Define rules to be run locally when using a compute cluster
localrules:
    _clairs_input_bam,
    _clairs_input_normal_bam,
    _clairs_get_resources,
    _clairs_link_clairs_models,
    _clairs_output_vcf,
    _clairs_all


VERSION_MAP_CLAIRS = CFG["options"]["version_map"]
SEQTYPE_MAP_CLAIRS = CFG["options"]["seqtype_map"]
MODULE_DIR = os.path.abspath(CFG["options"]["modsdir"])

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


def get_clairs_models(wc):
    CFG = config["lcr-modules"]["clairs"]
    base = CFG["options"]["model_path"]
    platform = get_platform(wc)
    models = {
        "full_alignment": f"{base}{platform}/full_alignment.pkl",
        "pileup": f"{base}{platform}/pileup.pkl"
        }
    if CFG["options"].get("call_indels", False):
        models["indel_full"] = f"{base}{platform}/indel/full_alignment.pkl"
        models["indel_pileup"] = f"{base}{platform}/indel/pileup.pkl"
    return models


# Clone ClairS repo at a fixed commit and download ClairS models
rule _clairs_get_resources:
    output:
        resources = touch(os.path.join(MODULE_DIR, "ClairS-0.4.4", "resources_dummy"))
    params:
        base_dir = os.path.join(MODULE_DIR, "ClairS-0.4.4"),
        repo_dir = os.path.join(MODULE_DIR, "ClairS-0.4.4", "ClairS"),
        models_dir = os.path.join(MODULE_DIR, "ClairS-0.4.4", "models"),
        tarball = os.path.join(MODULE_DIR, "ClairS-0.4.4", "models", "clairs_models.tar.gz"),
        commit = "4ff6c5fa3e59d5abb516ba4ed7d341d6b394197c"
    shell:
        op.as_one_line("""
            set -euo pipefail &&
            mkdir -p {params.base_dir} &&
            mkdir -p {params.models_dir} &&
            if [ ! -d {params.repo_dir}/.git ]; then git clone https://github.com/HKU-BAL/ClairS.git {params.repo_dir} ; fi &&
            cd {params.repo_dir} &&
            CURRENT_COMMIT="$(git rev-parse HEAD)" &&
            if [ "$CURRENT_COMMIT" != "{params.commit}" ]; then git fetch --all --tags && git checkout --force {params.commit} ; fi &&
            if ! find {params.models_dir} -type f -name 'pileup.pkl' | grep -q .; then wget -c -O {params.tarball} https://www.bio8.cs.hku.hk/clairs/models/clairs_models.tar.gz && tar -zxf {params.tarball} -C {params.models_dir} ; fi &&
            if ! find {params.models_dir} -type f -name 'pileup.pkl' | grep -q .; then echo "ERROR: No pileup.pkl found under models/" >&2 ; exit 1 ; fi &&
            if ! find {params.models_dir} -type f -name 'full_alignment.pkl' | grep -q .; then echo "ERROR: No full_alignment.pkl found under models/" >&2 ; exit 1 ; fi &&
            touch {output.resources}
        """)


# Link the ClairS models into the conda env bin
rule _clairs_link_clairs_models:
    input:
        resources = str(rules._clairs_get_resources.output.resources)
    output:
        models = touch(os.path.join(MODULE_DIR, "model_dummy"))
    log:
        stdout = CFG["logs"]["clairs"] + "link_clairs_models.stdout.log",
        stderr = CFG["logs"]["clairs"] + "link_clairs_models.stderr.log"
    conda:
        CFG["conda_envs"]["clairs"]
    params:
        link_target = lambda wc: os.path.join(config["lcr-modules"]["clairs"]["options"]["modsdir"], "ClairS-0.4.4", "models"),
        link_path = "$(dirname $(command -v pypy))/clairs_models"
    shell:
        op.as_one_line("""
            set -euo pipefail &&
            echo "link_target = {params.link_target}" >&2 &&
            echo "link_path   = {params.link_path}" >&2 &&
            TARGET="$(readlink -f {params.link_target})" &&
            if [ ! -d "$TARGET" ]; then echo "ERROR: Resolved ClairS models directory does not exist: $TARGET" >&2 ; exit 1 ; fi &&
            if ! find "$TARGET" -type f -name 'pileup.pkl' | grep -q .; then echo "ERROR: ClairS models directory exists but contains no pileup models" >&2 ; exit 1 ; fi &&
            if ! find "$TARGET" -type f -name 'full_alignment.pkl' | grep -q .; then echo "ERROR: ClairS models directory exists but contains no full-alignment models" >&2 ; exit 1 ; fi &&
            mkdir -p $(dirname {params.link_path}) &&
            ln -sfn "$TARGET" {params.link_path} &&
            test -d {params.link_path} &&
            touch {output.models}
        """)


# Calls variants using ClairS
rule _clairs_call_variants:
    input:
        tumour_bam = str(rules._clairs_input_bam.output.bam),
        normal_bam = str(rules._clairs_input_normal_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        models = str(rules._clairs_link_clairs_models.output.models)
    output:
        vcf = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.stdout.log",
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.stderr.log"
    params:
        clairs_args = CFG["options"]["clairs_args"],
        platform = get_platform,
        output_dir = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/",
        clairs_path = CFG["options"]["clairs_path"],
        models = lambda wc: get_clairs_models(wc),
        indel_p = lambda wc: (
            "--indel_pileup_model_path " + get_clairs_models(wc)["indel_pileup"]
            if "indel_pileup" in get_clairs_models(wc) else ""
            ),
        indel_fa = lambda wc: (
            "--indel_full_alignment_model_path " + get_clairs_models(wc)["indel_full"]
            if "indel_full" in get_clairs_models(wc) else ""
            ),
        call_indels = lambda wc: (
            "--enable_indel_calling --indel_min_qual 5"
            if config["lcr-modules"]["clairs"]["options"].get("call_indels", False) else ""
            )
    conda:
        CFG["conda_envs"]["clairs"]
    threads:
        CFG["threads"]["clairs"]
    resources:
        **CFG["resources"]["clairs"]
    shell:
        op.as_one_line("""
        {params.clairs_path}
            -P {params.models[pileup]}
            -F {params.models[full_alignment]}
            --tumor_bam_fn {input.tumour_bam}
            --normal_bam_fn {input.normal_bam}
            --ref_fn {input.fasta}
            --threads {threads}
            --platform {params.platform}
            --output_dir {params.output_dir}
            {params.call_indels}
            {params.indel_p}
            {params.indel_fa}
            {params.clairs_args}
            > {log.stdout} 2> {log.stderr}
        """)


# Combines ClairS VCF output files so indels and SNVs are in the same file
rule _clairs_combine_vcfs:
    input:
        vcf = str(rules._clairs_call_variants.output.vcf),
        tbi = str(rules._clairs_call_variants.output.tbi)
    output:
        combined_vcf = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combined.vcf.gz"),
        combined_tbi = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combined.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combine_vcfs.stdout.log",
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/combine_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
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
        bcftools sort {params.output_dir}/output.vcf.gz -Oz -o {params.output_dir}/output.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf {params.output_dir}/output.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        bcftools sort {params.output_dir}/indel.vcf.gz -Oz -o {params.output_dir}/indel.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf  {params.output_dir}/indel.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        bcftools concat
            --allow-overlaps
            {params.output_dir}/output.sorted.vcf.gz
            {params.output_dir}/indel.sorted.vcf.gz
            -Oz -o {output.combined_vcf} >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf {output.combined_vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Annotates VCF file with gnomAD frequency data and filters out poor calls
rule _clairs_gnomad_annotation:
    input:
        vcf = str(rules._clairs_combine_vcfs.output.combined_vcf),
        tbi = str(rules._clairs_combine_vcfs.output.combined_tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.gnomad.vcf.gz"),
        tbi = temp(CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/output.gnomad.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs_gnomad_annotation.stderr.log",
        stdout = CFG["logs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs_gnomad_annotation.stdout.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} 
        -a {input.gnomad} -c INFO/AF {input.vcf} 2> {log.stderr} |
        awk 'BEGIN {{FS=OFS="\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' |
        bcftools view -Oz -o {output.vcf} 2> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Filters out PoN variants
rule _clairs_filter:
    input:
        vcf = str(rules._clairs_gnomad_annotation.output.vcf), 
        tbi = str(rules._clairs_gnomad_annotation.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.final.vcf.gz",
        tbi = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/clairs.final.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/filter_pon.stderr.log",
        stdout = CFG["logs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/filter_pon.stdout.log"
    params:
        filters = CFG["options"]["filters"]
    shell:
        op.as_one_line("""
        bcftools isec -C -w1 {input.vcf} {input.pon} 2> {log.stderr} | 
        bcftools view -i '{params.filters}' -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _clairs_output_vcf:
    input:
        vcf = str(rules._clairs_filter.output.vcf),
        tbi = str(rules._clairs_filter.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched.clairs.combined.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched.clairs.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Cleans up additional files created by ClairS
rule _clairs_clean:
    input:
        str(rules._clairs_output_vcf.output.vcf)
    output:
        cleanup_complete = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}--unmatched/cleanup_complete.txt"
    shell:
        op.as_one_line("""
        d=$(dirname {output.cleanup_complete});
        rm -rf $d/tmp/ &&
        rm -rf $d/logs/ &&
        rm -f $d/run_clairs.log* &&
        rm -f $d/snv.vcf* &&
        rm -f $d/indel.* &&
        rm -f $d/output.* &&
        touch {output.cleanup_complete}
        """)


# Generates the target sentinels for each run, which generate the symlinks
rule _clairs_all:
    input:
        expand(
            [
                str(rules._clairs_output_vcf.output.vcf),
                str(rules._clairs_output_vcf.output.tbi),
                str(rules._clairs_clean.output.cleanup_complete)
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
