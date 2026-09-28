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
# `CFG` is a shortcut to `config["lcr-modules"]["clairs_to"]`
CFG = op.setup_module(
    name = "clairs_to",
    version = "1.1",
    subdirectories = ["inputs", "clairs_to", "filter", "gnomad", "outputs"]
)

config["pipeline_name"] = "clairs_to.yaml"

# Define rules to be run locally when using a compute cluster
localrules:
    _clairs_to_input_bam,
    _clairs_to_input_chrs,
    _clairs_to_chrom_bed,
    _clairs_to_clean,
    _clairs_to_cleanup_all_chroms,
    _clairs_to_output_vcf,
    _clairs_to_all


VERSION_MAP_CLAIRS_TO = CFG["options"]["version_map"]
SEQTYPE_MAP_CLAIRS_TO = CFG["options"]["seqtype_map"]

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


# Symlink chromosomes used for parallelization. genome_build-keyed only (not per-sample)
# Override via options.chromosomes_file (default excludes chrY)
CLAIRS_TO_CHROMOSOMES_FILE = CFG["options"].get("chromosomes_file", "")

checkpoint _clairs_to_input_chrs:
    input:
        chrs = CLAIRS_TO_CHROMOSOMES_FILE if CLAIRS_TO_CHROMOSOMES_FILE else reference_files("genomes/{genome_build}/genome_fasta/main_chromosomes.txt")
    output:
        chrs = CFG["dirs"]["inputs"] + "chroms/{genome_build}/main_chromosomes.txt"
    run:
        op.absolute_symlink(input.chrs, output.chrs)


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


# Optional restriction to a target-regions BED
# --bed_fn is ClairS-TO's only region-restriction flag, so it's reused for both the target panel and the chromosome split; no target BED means whole-chromosome, built from genome.fa.fai
CLAIRS_TO_TARGET_BED = CFG["options"].get("target_regions_bed", "")


rule _clairs_to_chrom_bed:
    input:
        bed = CLAIRS_TO_TARGET_BED if CLAIRS_TO_TARGET_BED else reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        bed = CFG["dirs"]["inputs"] + "chroms/{genome_build}/target_regions/{chrom}.bed"
    params:
        # fai branch synthesizes start/end (0, chrom_length) since a .fai row is name+length, not already BED-shaped
        using_fai = not bool(CLAIRS_TO_TARGET_BED)
    shell:
        op.as_one_line("""
        if [ "{params.using_fai}" = "True" ]; then
            awk -v chrom="{wildcards.chrom}" 'BEGIN {{FS=OFS="\\t"}} $1 == chrom {{print $1, 0, $2}}' {input.bed} > {output.bed}
        else
            awk -v chrom="{wildcards.chrom}" '$1 == chrom' {input.bed} > {output.bed}
        fi
        """)


def _clairs_to_get_chrom_bed(wildcards):
    CFG = config["lcr-modules"]["clairs_to"]
    return CFG["dirs"]["inputs"] + f"chroms/{wildcards.genome_build}/target_regions/{wildcards.chrom}.bed"


# Base dir for run_clairs_to's own --output_dir (raw per-chromosome scratch: tmp/, tmp_TUMOR/, logs/, snv_TUMOR.vcf.gz, indel_TUMOR.vcf.gz before combining)
# Defaults to the module's own clairs_to/ subdirectory; override via options.intermediate_results_dir_base to redirect onto different storage. Must be shared/network-visible, NOT node-local; same reasoning as deepsomatic/1.0's equivalent option. _clairs_to_combine_vcfs's own final combined.vcf.gz always stays module-managed regardless, so downstream rules are unaffected by this override
CLAIRS_TO_INTERMEDIATE_BASE = CFG["options"].get("intermediate_results_dir_base", "") or CFG["dirs"]["clairs_to"]


# Calls variants using ClairS-TO, one chromosome at a time via --bed_fn. output_dir is keyed by {chrom} too, so concurrent per-chromosome jobs don't clobber each other's files
rule _clairs_to_call_variants:
    input:
        tumour_bam = str(rules._clairs_to_input_bam.output.bam),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        chrom_bed = _clairs_to_get_chrom_bed
    output:
        indel_vcf = temp(CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/indel_TUMOR.vcf.gz"),
        indel_tbi = temp(CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/indel_TUMOR.vcf.gz.tbi"),
        snv_vcf = temp(CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/snv_TUMOR.vcf.gz"),
        snv_tbi = temp(CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/snv_TUMOR.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/{chrom}.clairs_to.stdout.log",
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/{chrom}.clairs_to.stderr.log"
    params:
        clairs_to_args = CFG["options"]["clairs_to_args"],
        platform = get_platform,
        output_dir = CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/"
    container:
        CFG["container_envs"]["clairs_to"]
    threads:
        CFG["threads"]["clairs_to_run"]
    resources:
        **CFG["resources"]["clairs_to_run"]
    shell:
        op.as_one_line("""
        /opt/bin/run_clairs_to
            --tumor_bam_fn {input.tumour_bam}
            --ref_fn {input.fasta}
            --bed_fn {input.chrom_bed}
            --threads {threads}
            --platform {params.platform}
            --output_dir {params.output_dir}
            --conda_prefix /opt/micromamba/envs/clairs-to
            -s TUMOR
            {params.clairs_to_args}
            >> {log.stdout} 2>> {log.stderr}
        """)


# Combines one chromosome's ClairS-TO VCF output files so indels and SNVs are in the same file
rule _clairs_to_combine_vcfs:
    input:
        indel_vcf = str(rules._clairs_to_call_variants.output.indel_vcf),
        indel_tbi = str(rules._clairs_to_call_variants.output.indel_tbi),
        snv_vcf = str(rules._clairs_to_call_variants.output.snv_vcf),
        snv_tbi = str(rules._clairs_to_call_variants.output.snv_tbi)
    output:
        vcf = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/combined.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/combined.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/{chrom}.combine_vcfs.stdout.log",
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/{chrom}.combine_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    params:
        output_dir = CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/"
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


# Cleans up one chromosome's additional ClairS-TO files as soon as that chromosome's own combined VCF has been written
rule _clairs_to_clean:
    input:
        str(rules._clairs_to_combine_vcfs.output.vcf)
    output:
        cleanup_complete = touch(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/cleanup_complete.txt")
    params:
        cleanup = CFG["options"]["cleanup_toggle"],
        # explicit path, not derived from output.cleanup_complete's own dirname; that's always module-managed, but the scratch files being cleaned up here live under CLAIRS_TO_INTERMEDIATE_BASE, which may be a different, user-redirected location
        scratch_dir = CLAIRS_TO_INTERMEDIATE_BASE + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/chromosomes/{chrom}/"
    shell:
        op.as_one_line("""
        d="{params.scratch_dir}" &&
        if [ "{params.cleanup}" = "True" ] || [ "{params.cleanup}" = "true" ]; then
            rm -rf "$d/tmp_TUMOR/" &&
            rm -rf "$d/tmp/" &&
            rm -rf "$d/logs/" &&
            rm -f "$d"/run_clairs_to.log* &&
            rm -f "$d"/snv.sorted.vcf.gz* &&
            rm -f "$d"/indel.sorted.vcf.gz* ;
        else
            echo "cleanup_toggle is false; Skipping cleanup" >&2 ;
        fi
        """)


def _clairs_to_get_chr_vcfs(wildcards):
    CFG = config["lcr-modules"]["clairs_to"]
    chrs = checkpoints._clairs_to_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    return expand(
        CFG["dirs"]["clairs_to"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{chemistry}}--tumor_only/chromosomes/{chrom}/combined.vcf.gz",
        chrom = chrs
    )


def _clairs_to_get_chr_cleanup(wildcards):
    CFG = config["lcr-modules"]["clairs_to"]
    chrs = checkpoints._clairs_to_input_chrs.get(**wildcards).output.chrs
    with open(chrs) as file:
        chrs = file.read().rstrip("\n").split("\n")
    return expand(
        CFG["dirs"]["clairs_to"] + "{{seq_type}}--{{genome_build}}/{{tumour_id}}--{{chemistry}}--tumor_only/chromosomes/{chrom}/cleanup_complete.txt",
        chrom = chrs
    )


# Merge per-chromosome combined VCFs into one sorted, indexed VCF at the same final path _clairs_to_combine_vcfs used to produce directly
rule _clairs_to_merge_vcfs:
    input:
        vcf = _clairs_to_get_chr_vcfs,
        cleanup = _clairs_to_get_chr_cleanup
    output:
        vcf = CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combined.vcf.gz",
        tbi = CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/combined.vcf.gz.tbi"
    log:
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to_merge_vcfs.stdout.log",
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to_merge_vcfs.stderr.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    threads:
        CFG["threads"]["clairs_to_merge_vcfs"]
    resources:
        **CFG["resources"]["clairs_to_merge_vcfs"]
    params:
        mem_mb = lambda wildcards, resources: int(resources.mem_mb * 0.8)
    shell:
        op.as_one_line("""
        bcftools concat --threads {threads} -a -O z {input.vcf} 2> {log.stderr}
            |
        bcftools sort -m {params.mem_mb}M -O z -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} > {log.stdout} 2>> {log.stderr}
        """)


# Aggregates all per-chromosome cleanups for one sample-run into a single dummy target
rule _clairs_to_cleanup_all_chroms:
    input:
        _clairs_to_get_chr_cleanup
    output:
        cleanup_complete = touch(CFG["dirs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/cleanup_complete.txt")


# Filters out poor quality variants
rule _clairs_to_filter:
    input:
        vcf = str(rules._clairs_to_merge_vcfs.output.vcf),
        tbi = str(rules._clairs_to_merge_vcfs.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.final.vcf.gz",
        tbi = CFG["dirs"]["filter"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_to.final.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    container:
        CFG["container_envs"]["bcftools"]
    resources:
        **CFG["resources"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/filter.stderr.log",
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/filter.stdout.log"
    params:
        min_depth = CFG["options"]["filters"]["min_depth"],
        snv_min_af = CFG["options"]["filters"]["snv"]["min_af"],
        snv_min_alt_depth = CFG["options"]["filters"]["snv"]["min_alt_depth"],
        indel_min_af = CFG["options"]["filters"]["indel"]["min_af"],
        indel_min_alt_depth = CFG["options"]["filters"]["indel"]["min_alt_depth"]
    shell:
        op.as_one_line("""
        bcftools isec -C -w1 {input.vcf} {input.pon} 2> {log.stderr} |
        bcftools view
            -i 'FILTER="PASS" &&
                FMT/DP[0] >= {params.min_depth} &&
                (
                    (TYPE="snp" && FMT/AF[0] >= {params.snv_min_af} && FMT/AD[0:1] >= {params.snv_min_alt_depth}) ||
                    (TYPE="indel" && FMT/AF[0] >= {params.indel_min_af} && FMT/AD[0:1] >= {params.indel_min_alt_depth})
                )'
            -Oz -o {output.vcf} 2>> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Annotates VCF file with gnomAD frequency data and filters out poor calls
rule _clairs_to_gnomad_annotation:
    input:
        vcf = str(rules._clairs_to_filter.output.vcf),
        tbi = str(rules._clairs_to_filter.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/output.gnomad.vcf.gz",
        tbi = CFG["dirs"]["gnomad"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/output.gnomad.vcf.gz.tbi"
    log:
        stderr = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_gnomad_annotation.stderr.log",
        stdout = CFG["logs"]["clairs_to"] + "{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only/clairs_gnomad_annotation.stdout.log"
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
rule _clairs_to_output_vcf:
    input:
        vcf = str(rules._clairs_to_gnomad_annotation.output.vcf),
        tbi = str(rules._clairs_to_gnomad_annotation.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only.clairs_to.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{chemistry}--tumor_only.clairs_to.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _clairs_to_all:
    input:
        expand(
            [
                str(rules._clairs_to_output_vcf.output.vcf),
                str(rules._clairs_to_output_vcf.output.tbi),
                str(rules._clairs_to_cleanup_all_chroms.output.cleanup_complete)
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
