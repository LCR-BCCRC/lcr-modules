#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Giuliano Banco
# Module Author:    Giuliano Banco
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
# `CFG` is a shortcut to `config["lcr-modules"]["clairs_to"]`
CFG = op.setup_module(
    name = "clairs_to",
    version = "1.0",
    subdirectories = ["inputs", "clairs_to", "clairs_to_gnomad", "clairs_to_filter", "outputs"],
)

config["pipeline_name"] = "clairs_to.yaml"

# Define rules to be run locally when using a compute cluster
localrules:
    _clairs_to_input_bam,
    _clairs_to_get_resources,
    _clairs_to_link_clairs_to_models,
    _clairs_to_output_vcf,
    _clairs_to_all


VERSION_MAP_CLAIRS_TO = CFG["options"]["version_map"]
SEQTYPE_MAP_CLAIRS_TO = CFG["options"]["seqtype_map"]
MODULE_DIR = os.path.abspath(CFG["options"]["modsdir"])

possible_genome_builds = ", ".join(list(VERSION_MAP_CLAIRS_TO.keys()))
for genome_build in CFG["runs"]["tumour_genome_build"]:
    assert genome_build in possible_genome_builds, (
        f"Samples table includes genome builds not yet compatible with this module. "
        f"This module is currently only compatible with {possible_genome_builds}. "
    )

possible_seq_types = ", ".join(list(SEQTYPE_MAP_CLAIRS_TO.keys()))
for seq_type in CFG["runs"]["tumour_seq_type"]:
    assert seq_type in possible_seq_types, (
        f"Samples table includes seq types not yet compatible with this module. "
        f"This module is currently only compatible with {possible_seq_types}. "
    )


chemistry = CFG["runs"]["tumour_chemistry"]


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _clairs_to_input_bam:
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


def get_platform(wildcards):
    CFG = config["lcr-modules"]["clairs_to"]
    this_sample = op.filter_samples(
        CFG["runs"],
        sample_id = wildcards.tumour_id,
        seq_type = wildcards.seq_type, 
        genome_build = wildcards.genome_build
    )
    platform = this_sample["tumour_platform"].tolist()[0]
    return platform


def get_clairs_to_models(wc):
    CFG = config["lcr-modules"]["clairs_to"]
    base = CFG["options"]["model_path"]
    platform = get_platform(wc)
    models = {
        "snv_affirmative": f"{base}{platform}/pileup_affirmative.pkl",
        "snv_negational": f"{base}{platform}/pileup_negational.pkl",
        "indel_affirmative": f"{base}{platform}/indel/pileup_affirmative.pkl",
        "indel_negational": f"{base}{platform}/indel/pileup_negational.pkl"
        }
    return models


# Clone ClairS-TO repo at a fixed commit and download ClairS-TO resources
rule _clairs_to_get_resources:
    output:
        resources = touch(os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources_dummy"))
    params:
        base_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2"),
        repo_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "ClairS-TO"),
        resources_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources"),
        models_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_models"),
        databases_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_databases"),
        references_dir = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_cna_data"),
        models_tarball = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_models", "clairs-to_models.tar.gz"),
        databases_tarball = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_databases", "clairs-to_databases.tar.gz"),
        references_tarball = os.path.join(MODULE_DIR, "ClairS-TO-0.4.2", "resources", "clairs-to_cna_data", "clairs-to_reference_files.tar.gz"),
        commit = "1179868508ca47d7972fc5d2b3fb644c0e66845c"
    shell:
        op.as_one_line("""
        set -euo pipefail
        &&
        mkdir -p {params.base_dir} {params.resources_dir} {params.models_dir} {params.databases_dir} {params.references_dir} &&
        if [ ! -d {params.repo_dir}/.git ]; then git clone https://github.com/HKU-BAL/ClairS-TO.git {params.repo_dir} ; fi &&
        cd {params.repo_dir} &&
        CURRENT_COMMIT="$(git rev-parse HEAD)" &&
        if [ "$CURRENT_COMMIT" != "{params.commit}" ]; then git fetch --all --tags && git checkout --force {params.commit} ; fi &&
        if ! find {params.models_dir} -type f -name 'pileup_affirmative.pkl' -print -quit | read ; then
            wget -c -O {params.models_tarball} http://www.bio8.cs.hku.hk/clairs-to/models/clairs-to_models.tar.gz &&
            tar -zxf {params.models_tarball} -C {params.models_dir} ; fi &&
        if ! find {params.databases_dir} -type f -name '1000g-pon.sites.vcf.gz' -print -quit | read ; then
            wget -c -O {params.databases_tarball} http://www.bio8.cs.hku.hk/clairs-to/databases/clairs-to_databases.tar.gz &&
            tar -zxf {params.databases_tarball} -C {params.databases_dir} ; fi &&
        if ! find {params.references_dir} -type f -name 'GC_G1000_hg38.txt' -print -quit | read ; then
            wget -c -O {params.references_tarball} http://www.bio8.cs.hku.hk/clairs-to/cna_data/reference_files.tar.gz &&
            tar -zxf {params.references_tarball} -C {params.references_dir} ;
        fi &&
        if ! find {params.models_dir} -type f -name 'pileup_affirmative.pkl' -print -quit | read ; then
            echo "ERROR: Missing pileup_affirmative.pkl in models directory=" >&2 ; exit 1 ;
        fi
        &&
        touch {output.resources}
        """)


# Link the ClairS-TO models into the conda env bin
rule _clairs_to_link_clairs_to_models:
    input:
        resources = str(rules._clairs_to_get_resources.output.resources)
    output:
        link_dummy = touch(os.path.join(MODULE_DIR, "link_dummy"))
    conda:
        CFG["conda_envs"]["clairs_to"]
    params:
        models_target = lambda wc: os.path.join(config["lcr-modules"]["clairs_to"]["options"]["modsdir"], "ClairS-TO-0.4.2", "resources", "clairs-to_models"),
        models_link   = "$(dirname $(command -v python))/clairs-to_models",
        databases_target = lambda wc: os.path.join(config["lcr-modules"]["clairs_to"]["options"]["modsdir"], "ClairS-TO-0.4.2", "resources", "clairs-to_databases"),
        databases_link   = "$(dirname $(command -v python))/clairs-to_databases",
        references_target = lambda wc: os.path.join(config["lcr-modules"]["clairs_to"]["options"]["modsdir"], "ClairS-TO-0.4.2", "resources", "clairs-to_cna_data"),
        references_link   = "$(dirname $(command -v python))/clairs-to_cna_data"
    shell:
        op.as_one_line("""
        set -euo pipefail
        &&
        TARGET=$(readlink -f {params.models_target}) &&
        LINK={params.models_link} &&
        mkdir -p $(dirname "$LINK") &&
        ln -sfn "$TARGET" "$LINK" &&
        test -d "$LINK"
        &&
        TARGET=$(readlink -f {params.databases_target}) &&
        LINK={params.databases_link} &&
        mkdir -p $(dirname "$LINK") &&
        ln -sfn "$TARGET" "$LINK" &&
        test -d "$LINK"
        &&
        TARGET=$(readlink -f {params.references_target}) &&
        LINK={params.references_link} &&
        mkdir -p $(dirname "$LINK") &&
        ln -sfn "$TARGET" "$LINK" &&
        test -d "$LINK"
        &&
        touch {output.link_dummy}
        """)


# Calls variants using ClairS-TO
rule _clairs_to_call_variants:
    input:
        tumour_bam = str(rules._clairs_to_input_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        models = str(rules._clairs_to_link_clairs_to_models.output.link_dummy)
    output:
        indel_vcf = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/indel.vcf.gz"),
        indel_tbi = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/indel.vcf.gz.tbi"),
        snv_vcf = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/snv.vcf.gz"),
        snv_tbi = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/snv.vcf.gz.tbi")        
    log:
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.stdout.log",
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.stderr.log"
    params:
        clairs_to_args = CFG["options"]["clairs_to_args"],
        platform = get_platform,
        output_dir = CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/",
        clairs_to_path = CFG["options"]["clairs_to_path"],
        models = lambda wc: get_clairs_to_models(wc)
    conda:
        CFG["conda_envs"]["clairs_to"]
    threads:
        CFG["threads"]["clairs_to"]
    resources:
        **CFG["resources"]["clairs_to"]
    shell:
        op.as_one_line("""
        {params.clairs_to_path}
            --snv_pileup_affirmative_model_path {params.models[snv_affirmative]}
            --snv_pileup_negational_model_path {params.models[snv_negational]}
            --indel_pileup_affirmative_model_path {params.models[indel_affirmative]}
            --indel_pileup_negational_model_path {params.models[indel_negational]}
            --tumor_bam_fn {input.tumour_bam}
            --ref_fn {input.fasta}
            --threads {threads}
            --platform {params.platform}
            --output_dir {params.output_dir}
            {params.clairs_to_args}
            > {log.stdout} 2> {log.stderr}
        """)


# Combines ClairS VCF output files so indels and SNVs are in the same file
rule _clairs_to_combine_vcfs:
    input:
        indel_vcf = str(rules._clairs_to_call_variants.output.indel_vcf),
        indel_tbi = str(rules._clairs_to_call_variants.output.indel_tbi),
        snv_vcf = str(rules._clairs_to_call_variants.output.snv_vcf),
        snv_tbi = str(rules._clairs_to_call_variants.output.snv_tbi)
    output:
        vcf = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combined.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combined.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combine_vcfs.stdout.log",
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combine_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    params:
        output_dir = CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/",
    shell:
        op.as_one_line("""
        bcftools sort {input.snv_vcf} -Oz -o {params.output_dir}/snv.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf {params.output_dir}/snv.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        bcftools sort {input.indel_vcf} -Oz -o {params.output_dir}/indel.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf  {params.output_dir}/indel.sorted.vcf.gz >> {log.stdout} 2>> {log.stderr} &&
        bcftools concat
            --allow-overlaps
            {params.output_dir}/snv.sorted.vcf.gz
            {params.output_dir}/indel.sorted.vcf.gz
            -Oz -o {output.vcf} >> {log.stdout} 2>> {log.stderr} &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Annotates VCF file with gnomAD frequency data and filters out poor calls
rule _clairs_to_gnomad_annotation:
    input:
        vcf = str(rules._clairs_to_combine_vcfs.output.vcf),
        tbi = str(rules._clairs_to_combine_vcfs.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["clairs_to_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/output.gnomad.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs_to_gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/output.gnomad.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_gnomad_annotation.stderr.log",
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_gnomad_annotation.stdout.log"
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
rule _clairs_to_filter:
    input:
        vcf = str(rules._clairs_to_gnomad_annotation.output.vcf), 
        tbi = str(rules._clairs_to_gnomad_annotation.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["clairs_to_filter"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.final.vcf.gz",
        tbi = CFG["dirs"]["clairs_to_filter"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.final.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/filter_pon.stderr.log",
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/filter_pon.stdout.log"
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
rule _clairs_to_output_vcf:
    input:
        vcf = str(rules._clairs_to_filter.output.vcf),
        tbi = str(rules._clairs_to_filter.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only.clairs_to.combined.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only.clairs_to.combined.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Cleans up additional files created by ClairS
rule _clairs_to_clean:
    input:
        str(rules._clairs_to_output_vcf.output.vcf)
    output:
        cleanup_complete = CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/cleanup_complete.txt"
    shell:
        op.as_one_line("""
        d=$(dirname {output.cleanup_complete});
        rm -rf $d/tmp/ &&
        rm -rf $d/logs/ &&
        rm -f $d/run_clairs_to.log* &&
        rm -f $d/snv.* &&
        rm -f $d/indel.* &&
        touch {output.cleanup_complete}
        """)


# Generates the target sentinels for each run, which generate the symlinks
rule _clairs_to_all:
    input:
        expand(
            [
                str(rules._clairs_to_output_vcf.output.vcf),
                str(rules._clairs_to_output_vcf.output.tbi),
                str(rules._clairs_to_clean.output.cleanup_complete)
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            platform=CFG["runs"]["tumour_platform"],
            chemistry=CFG["runs"]["tumour_chemistry"]
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
