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
    subdirectories = ["inputs", "clairs", "outputs"],
)

config["pipeline_name"] = "clairs.yaml"

# Define rules to be run locally when using a compute cluster
localrules:
    _clairs_input_bam,
    _clairs_output_vcf,
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

if "tumour_chemistry" not in CFG["runs"].columns:
    raise KeyError("Chemistry is needed in your samples table. Available chemistry values: R9, R10.")

if "tumour_platform" not in CFG["runs"].columns:
    raise KeyError("Platform is needed in your samples table. See config for possible values.")

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

# Rule to see if current container version is most recent?

rule _clairs_install:
    output:
        sif = CFG["dirs"]["clairs"] + "clairs.sif"
    log:
        stdout = CFG["logs"]["clairs"] + "install.stdout.log",
        stderr = CFG["logs"]["clairs"] + "install.stderr.log"
    conda:
        CFG["conda_envs"]["singularity"]
    shell:
        """
        singularity pull {output.sif} docker://hkubal/clairs:latest
        """


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
        tumour_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        normal_bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{normal_name}.bam",
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        fai = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai"),
        sif = str(rules._clairs_install.output.sif)
    output:
        vcf = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/output.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/output.vcf.gz.tbi")
    log:
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs.stdout.log",
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs.stderr.log"
    params:
        options = CFG["options"]["clairs"],
        container = CFG["dirs"]["clairs"] + "clairs.sif",
        platform = get_platform,
        output_dir = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/",
        tumor_bam_real = lambda wc, input: os.path.realpath(input.tumour_bam),
        normal_bam_real = lambda wc, input: os.path.realpath(input.normal_bam),
        ref_fasta_real = lambda wc, input: os.path.realpath(input.fasta),
        tumor_dir = lambda wc, input: os.path.dirname(os.path.realpath(input.tumour_bam)),
        normal_dir = lambda wc, input: os.path.dirname(os.path.realpath(input.normal_bam)),
        ref_dir = lambda wc, input: os.path.dirname(os.path.realpath(input.fasta)),
        tumor_base = lambda wc, input: os.path.basename(os.path.realpath(input.tumour_bam)),
        normal_base = lambda wc, input: os.path.basename(os.path.realpath(input.normal_bam)),
        ref_base = lambda wc, input: os.path.basename(os.path.realpath(input.fasta))
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["clairs"]
    resources:
        **CFG["resources"]["clairs"]
    shell:
        op.as_one_line("""
            apptainer exec
                -B {params.tumor_dir}:/tumor:ro
                -B {params.normal_dir}:/normal:ro
                -B {params.ref_dir}:/ref:ro
                -B {params.output_dir}:/output
                {params.container} bash -c "
                    /opt/bin/run_clairs
                        --tumor_bam_fn /tumor/{params.tumor_base}
                        --normal_bam_fn /normal/{params.normal_base}
                        --ref_fn /ref/{params.ref_base}
                        --threads {threads}
                        --output_dir /output
                        -p {params.platform}
                        --conda_prefix /opt/conda/envs/clairs
                        {params.options}
                "
                > {log.stdout} 2> {log.stderr} 
        """)


    # Annotates VCF file with gnomAD frequency data and filters out poor calls.
rule _clairs_gnomad_annotation:
    input:
        vcf = str(rules._clairs_call_variants.output.vcf),
        tbi = str(rules._clairs_call_variants.output.tbi),
        gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/output.gnomad.vcf.gz"),
        tbi = temp(CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/output.gnomad.vcf.gz.tbi")
    log:
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs_gnomad_annotation.stderr.log",
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs_gnomad_annotation.stdout.log"
    conda:
        CFG["conda_envs"]["bcftools"]
    threads:
        CFG["threads"]["bcftools"]
    resources:
        **CFG["resources"]["gnomad"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads}
        -a {input.gnomad} -c INFO/AF {input.vcf} |
        awk 'BEGIN {{FS=OFS="\\t"}} {{ if ($1 !~ /^#/ && $8 !~ ";AF=") $8=$8";AF=0"; print $0; }}' |
        bcftools view -i 'FILTER="PASS" && INFO/AF < 0.0001 && FMT/DP[0] >= 8 && FMT/AF[0] >= 0.1 && (FMT/NDP[0] < 10 || FMT/NAF[0] < 0.3)'
        -Oz -o {output.vcf} 2> {log.stderr}
        &&
        tabix -p vcf {output.vcf} 2>> {log.stderr}
        """)


# Filters out PoN variants
rule _clairs_filter_pon:
    input:
        vcf = str(rules._clairs_gnomad_annotation.output.vcf), 
        tbi = str(rules._clairs_gnomad_annotation.output.tbi),
        pon = reference_files("genomes/{genome_build}/ont/colorsDb.v1.2.0.deepvariant.glnexus.{genome_build}.vcf.gz")
    output:
        vcf = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs.final.vcf.gz",
        tbi = CFG["dirs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/clairs.final.vcf.gz.tbi"
    conda:
        CFG["conda_envs"]["bcftools"]
    resources:
        **CFG["resources"]["filter"]
    threads:
        CFG["threads"]["bcftools"]
    log:
        stderr = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/filter_pon.stderr.log",
        stdout = CFG["logs"]["clairs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}/filter_pon.stdout.log"
    shell:
        op.as_one_line("""
        bcftools isec -C -w1 -O z -o {output.vcf} {input.vcf} {input.pon} > {log.stdout} 2> {log.stderr}
        &&
        tabix -p vcf {output.vcf} >> {log.stdout} 2>> {log.stderr}
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _clairs_output_vcf:
    input:
        vcf = str(rules._clairs_filter_pon.output.vcf),
        tbi = str(rules._clairs_filter_pon.output.tbi)
    output:
        vcf = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.clairs.vcf.gz",
        tbi = CFG["dirs"]["outputs"] + "vcf/{seq_type}--{genome_build}/{tumour_id}--{normal_name}--{chemistry}.clairs.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf, output.vcf, in_module= True)
        op.relative_symlink(input.tbi, output.tbi, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _clairs_all:
    input:
        expand(
            [
                str(rules._clairs_output_vcf.output.vcf),
                str(rules._clairs_output_vcf.output.tbi)
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
