#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


##### SETUP #####

import sys, os
from os.path import join

import oncopipe as op
import pandas as pd

# Setup module and store module-specific configuration in `CONFIG`
CFG = op.setup_module(
    name = "vcf2maf",
    version = "1.3",
    subdirectories = ["inputs","decompressed","vcf2maf","crossmap", "normalize", "outputs"]
)

# Define rules to be run locally when using a compute cluster
localrules:
    _vcf2maf_input_vcf,
    _vcf2maf_gnomad_filter_maf,
    _vcf2maf_output_maf,
    _vcf2maf_crossmap,
    _vcf2maf_all

#set variable for prepending to PATH based on config
VCF2MAF_SCRIPT_PATH = CFG['inputs']['src_dir']

# This is going to be annoying, but to normalize prefixes, we need to determine
# 1) Which iteration of the genome build something is (i.e. GRch37 or GRCh38)
# 2) If a genome build is chr-prefixed or not
# To determine this, lets load the reference config and re-parse it
VCF2MAF_REFERENCE_CONFIG = CFG["options"]["reference_config"]  # Reference config path
configfile: VCF2MAF_REFERENCE_CONFIG
# Store all the attributes we will need
VCF2MAF_GENOME_VERSION = {}  # Will be a simple {"GRCh38-SFU": "grch38"} etc.
VCF2MAF_GENOME_PREFIX = {}  # Will be a simple hash of {"GRCh38-SFU": True} if chr-prefixed
VCF2MAF_VERSION_MAP = {}

for genome_build, attributes in config['genome_builds'].items():
    try:
        genome_version = attributes["version"]
    except KeyError as e:  
        # This wasn't included in the reference entry for this genome build
        # This should never happen, as the reference workflow checks for this,
        # but ¯\_(ツ)_/¯
        raise AttributeError(f"Unable to determine the \"version\" of genome {genome_version} in reference config {VCF2MAF_REFERENCE_CONFIG}") from e
    try:
        genome_provider = attributes["provider"]
    except KeyError as e:
        raise AttributeError(f"Unable to determine the \"provider\" of genome {genome_version} in reference config {VCF2MAF_REFERENCE_CONFIG}") from e

    VCF2MAF_GENOME_VERSION[genome_build] = genome_version  # What is the parent genome build?
    VCF2MAF_GENOME_PREFIX[genome_build] = True if genome_provider == "ucsc" else False  # Is this chr-prefixed?
    VCF2MAF_VERSION_MAP[genome_build] = genome_version.replace("grch", "GRCh")  # Genome build for vcf2maf

##### RULES #####

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _vcf2maf_input_vcf:
    input:
        vcf_gz = CFG["inputs"]["sample_vcf_gz"]
    output:
        vcf_gz = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz",
        index = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.vcf.gz.tbi"
    run:
        op.relative_symlink(input.vcf_gz, output.vcf_gz)
        op.relative_symlink(input.vcf_gz + ".tbi", output.index)

rule _vcf2maf_annotate_gnomad:
    input:
        vcf = str(rules._vcf2maf_input_vcf.output.vcf_gz),
        normalized_gnomad = reference_files("genomes/{genome_build}/variation/af-only-gnomad.normalized.{genome_build}.vcf.gz")
    output:
        vcf = temp(CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.annotated.vcf")
    conda:
        CFG["conda_envs"]["bcftools"]
    resources: 
        **CFG["resources"]["annotate"]
    threads: 
        CFG["threads"]["annotate"]
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"]
    shell:
        op.as_one_line("""
        bcftools annotate --threads {threads} -a {input.normalized_gnomad} {input.vcf} -c "INFO/gnomADg_AF:=INFO/AF" -o {output.vcf}
        """)

rule _vcf2maf_run:
    input:
        vcf = str(rules._vcf2maf_annotate_gnomad.output.vcf),
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa"),
        vep_cache = CFG["inputs"]["vep_cache"]
    output:
        maf = temp(CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf"),
        vep = temp(CFG["dirs"]["decompressed"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.annotated.vep.vcf")
    log:
        stdout = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}_vcf2maf.stderr.log",
    params:
        opts = CFG["options"]["vcf2maf"],
        build = lambda w: VCF2MAF_VERSION_MAP[w.genome_build],
        custom_enst = lambda w: config['lcr-modules']["vcf2maf"]["switches"]["custom_enst"][VCF2MAF_GENOME_VERSION[w.genome_build]]
    conda:
        CFG["conda_envs"]["vcf2maf"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        **CFG["resources"]["vcf2maf"]
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"]
    shell:
        op.as_one_line("""
        VCF2MAF_SCRIPT_PATH={VCF2MAF_SCRIPT_PATH};
        PATH=$VCF2MAF_SCRIPT_PATH:$PATH;
        VCF2MAF_SCRIPT="$VCF2MAF_SCRIPT_PATH/vcf2maf.pl";
        if [[ -e {output.maf} ]]; then rm -f {output.maf}; fi;
        if [[ -e {output.vep} ]]; then rm -f {output.vep}; fi;
        vepPATH=$(dirname $(which vep))/../share/variant-effect-predictor*;
        if [[ $(which vcf2maf.pl) =~ $VCF2MAF_SCRIPT ]]; then
            echo "using bundled patched script $VCF2MAF_SCRIPT";
            echo "Using $VCF2MAF_SCRIPT to run {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
            vcf2maf.pl
            --input-vcf {input.vcf}
            --output-maf {output.maf}
            --tumor-id {wildcards.tumour_id}
            --normal-id {wildcards.normal_id}
            --ref-fasta {input.fasta}
            --ncbi-build {params.build}
            --vep-data {input.vep_cache}
            --vep-path $vepPATH
            {params.opts}
            --custom-enst {params.custom_enst}
            --retain-info gnomADg_AF
            >> {log.stdout} 2>> {log.stderr};
        else echo "WARNING: PATH is not set properly, using $(which vcf2maf.pl) will result in error during execution. Please ensure $VCF2MAF_SCRIPT exists." > {log.stderr};fi &&
        touch {output.vep}
        """)

rule _vcf2maf_gnomad_filter_maf:
    input:
        maf = str(rules._vcf2maf_run.output.maf)
    output:
        maf = temp(CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.maf"),
        dropped_maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf.gz"
    params:
        opts = CFG["options"]["gnomAD_cutoff"],
        temp_file = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf"
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"]
    shell:
        op.as_one_line("""
        cat {input.maf} | perl -lane 'next if /^(!?#)/; my @cols = split /\t/; @AF_all =split/,/, $cols[114]; $skip=0; for(@AF_all){{$skip++ if $_ > {params.opts}}} if ($skip) {{print STDERR;}} else {{print;}};' > {output.maf} 2>{params.temp_file}
            &&
        gzip {params.temp_file}
            &&
        touch {output.dropped_maf}
        """)

def get_original_genome(wildcards):
    # Determine the original (i.e. input) reference genome for this sample
    # Since this module projects to various output genome builds, we need to parse the sample table for the starting build
    # To determine what we need to do
    sample_table = config['lcr-modules']["vcf2maf"]["samples"]
    sample_entry = sample_table.loc[(sample_table["sample_id"] == wildcards.tumour_id) & (sample_table["seq_type"] == wildcards.seq_type)]
    if len(sample_entry) == 0:
        raise AttributeError("Unable to locate a a sample with tumour_id:{wildcards.tumour_id}, normal_id:{wildcards.normal_id}, seq_type:{wildcards.seq_type} in the \'sample\' table")
    original_genome_build = sample_entry.iloc[0]["genome_build"]
    return original_genome_build

def get_chain(genome_build):
    # NOTE: This only currently supports hg38 and hg19. If you are using other genome builds, this will need to be handled
    genome_version = VCF2MAF_GENOME_VERSION[genome_build]
    if genome_version == "grch38":
        return reference_files("genomes/" + genome_build + "/chains/grch38/hg38ToHg19.over.chain")
    elif genome_version == "grch37":
        return reference_files("genomes/" + genome_build +"/chains/grch37/hg19ToHg38.over.chain")
    else:
        raise AttributeError(f"No supported CrossMap chain for {genome_version} within this module")


def crossmap_input(wildcards):
    original_genome_build = get_original_genome(wildcards)
    return {"maf": expand(rules._vcf2maf_gnomad_filter_maf.output.maf, **wildcards, genome_build = original_genome_build),
            "convert_coord": config['lcr-modules']["vcf2maf"]["inputs"]["convert_coord"],
            "chains": get_chain(original_genome_build)}

rule _vcf2maf_crossmap:
    input:
        unpack(crossmap_input)
    output:
        maf = temp(CFG["dirs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf")
    log:
        stdout = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stdout.log",
        stderr = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.crossmap.stderr.log"
    conda:
        CFG["conda_envs"]["crossmap"]
    threads:
        CFG["threads"]["vcf2maf"]
    resources:
        **CFG["resources"]["crossmap"]
    wildcard_constraints:
        target_build = "hg38|hg19",  # Crossmap only converts to chr-prefixed outputs, so these are what will be generated
        base_name = CFG["vcf_base_name"]
    shell:
        op.as_one_line("""
        {input.convert_coord}
        {input.maf}
        {input.chains}
        {output.maf}
        crossmap
        > {log.stdout} 2> {log.stderr}
        """)


def get_normalize_input(wildcards, genome_build_only = False):
    new_genome_build  = wildcards.target_build
    # Since snakemake only knows what the TARGET genome build is, we need to find the source
    original_genome_build = get_original_genome(wildcards)
    original_genome_version = VCF2MAF_GENOME_VERSION[original_genome_build]
    new_genome_version = VCF2MAF_GENOME_VERSION[new_genome_build]

    # If using this as an input function for snakemake
    if not genome_build_only:
        # Do we need to run CrossMap on this? Check the genome version
        if original_genome_version != new_genome_version:
            # Source does not match, get the converted MAF from CrossMap
            wildcards.target_build = "hg38" if new_genome_version == "grch38" else "hg19"  # Since CrossMap only outputs to these types
            return expand(rules._vcf2maf_crossmap.output.maf, **wildcards)
        else:
            # Source matches. CrossMap not necessary
            return expand(rules._vcf2maf_gnomad_filter_maf.output.maf, **wildcards, genome_build = original_genome_build)
    else:
        # Just return the original genome build
        if original_genome_version != new_genome_version:
            return "hg38" if new_genome_version == "grch38" else "hg19"
        else:
            return original_genome_version

# Add or remove chr prefix as necessary
rule _vcf2maf_normalize_prefix:
    input:
        maf = get_normalize_input
    output:
        maf = CFG["dirs"]["normalize"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.maf"
    params:
        dest_chr = lambda w: VCF2MAF_GENOME_PREFIX[w.target_build]
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"]), 
        base_name = CFG["vcf_base_name"]
    run:
        input.maf = input.maf[0]  # Because of the input function and expand(), the input is technically a list of length 1
        maf_open = pd.read_csv(input.maf, sep = "\t")
        # Update build column in the MAF
        maf_open["NCBI_Build"] = wildcards.target_build
        # To handle CrossMap weirdness, remove all chr-prefixes and add them back later
        maf_open["Chromosome"] = maf_open["Chromosome"].astype(str).str.replace('chr', '')
        if params.dest_chr:  # Will evaluate to True if the destination genome is chr-prefixed
            # Add chr prefix
            maf_open['Chromosome'] = 'chr' + maf_open['Chromosome'].astype(str)

        maf_open.to_csv(output.maf, sep="\t", index=False)

rule _vcf2maf_output_maf_projection:
    input:
        maf = str(rules._vcf2maf_normalize_prefix.output.maf),
    output:
        maf = CFG["dirs"]["outputs"] + "{seq_type}--{target_build}--projection/{tumour_id}--{normal_id}--{pair_status}_{base_name}.maf"
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"]), 
        base_name = CFG["vcf_base_name"]
    run:
        op.relative_symlink(input.maf, output.maf)

rule _vcf2maf_output_maf:
    input:
        maf = str(rules._vcf2maf_normalize_prefix.output.maf),
    output:
        maf = CFG["dirs"]["outputs"] + "{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}_{base_name}.{target_build}.maf"
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"]), 
        base_name = CFG["vcf_base_name"]
    run:
        op.relative_symlink(input.maf, output.maf)

# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(
            expand(str(rules._vcf2maf_output_maf.output.maf), zip,
                seq_type = CFG["runs"]["tumour_seq_type"],
                tumour_id = CFG["runs"]["tumour_sample_id"],
                normal_id = CFG["runs"]["normal_sample_id"],
                pair_status = CFG["runs"]["pair_status"],
                base_name = [CFG["vcf_base_name"]] * len(CFG["runs"]["tumour_sample_id"]),
                allow_missing = True),
            target_build = CFG["options"]["target_builds"])
        # Why are there two expand statements? Well we want every iteration of these MAFs for all the target genome builds
        # But that is the only wildcard we want to expand to all iterations
        # Hence the inner expand is using zip()

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
