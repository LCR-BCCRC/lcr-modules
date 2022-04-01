#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Bruno Grande
# Module Author:    Helena Winata
# Contributors:     Ryan Morin


##### SETUP #####

import sys, os
from os.path import join
import hashlib

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
    _vcf2maf_input_bam, 
    _vcf2maf_output_maf,
    _vcf2maf_gnomad_filter_maf, 
    _vcf2maf_install_GAMBLR,
    _vcf2maf_output_original, 
    _vcf2maf_normalize_prefix,
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

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later. 
AUGMENT_SSM = os.path.abspath(config["lcr-modules"]["vcf2maf"]["inputs"]["augment_ssm"])

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
    group: "vcf_and_annotate"
    run:
        op.absolute_symlink(input.vcf_gz, output.vcf_gz)
        op.absolute_symlink(input.vcf_gz + ".tbi", output.index)

rule _vcf2maf_input_bam: 
    input: 
        bam = CFG["inputs"]["sample_bam"], 
        bai = CFG["inputs"]["sample_bai"]
    output: 
        bam = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/tumour.bam", 
        bai = CFG["dirs"]["inputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/tumour.bam.bai"
    group: "bam_and_augment"
    run: 
        op.absolute_symlink(input.bam, output.bam)
        op.absolute_symlink(input.bai, output.bai)

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
    group: "vcf_and_annotate"
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
        custom_enst = lambda w: "--custom-enst " + str(config['lcr-modules']["vcf2maf"]["switches"]["custom_enst"][VCF2MAF_GENOME_VERSION[w.genome_build]]) if config['lcr-modules']["vcf2maf"]["switches"]["custom_enst"][VCF2MAF_GENOME_VERSION[w.genome_build]] != "" else ""
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
        if [[ $(which vcf2maf.pl) =~ $(ls $VCF2MAF_SCRIPT) ]]; then
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
            {params.custom_enst}
            --retain-info gnomADg_AF
            >> {log.stdout} 2>> {log.stderr};
        else echo "ERROR: PATH is not set properly, using $(which vcf2maf.pl) will result in error during execution. Please ensure $VCF2MAF_SCRIPT exists." > {log.stderr}; exit 1; fi &&
        touch {output.vep}
        """)

rule _vcf2maf_gnomad_filter_maf:
    input:
        maf = str(rules._vcf2maf_run.output.maf)
    output:
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.raw.maf",
        dropped_maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf.gz"
    params:
        opts = CFG["options"]["gnomAD_cutoff"],
        temp_file = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.dropped.maf"
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"], 
        maf = "maf"
    shell:
        op.as_one_line("""
        cat {input.maf} | perl -lane 'next if /^(!?#)/; my @cols = split /\t/; @AF_all =split/,/, $cols[114]; $skip=0; for(@AF_all){{$skip++ if $_ > {params.opts}}} if ($skip) {{print STDERR;}} else {{print;}};' > {output.maf} 2>{params.temp_file}
            &&
        gzip {params.temp_file}
            &&
        touch {output.dropped_maf}
        """)

# Obtain the path to the GAMBLR conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix: 
    conda_prefix = workflow.conda_prefix
else: 
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open("config/envs/GAMBLR.yaml", 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
GAMBLR = glob.glob(conda_prefix + "/" + h[:8] + "*")[0]

rule _vcf2maf_install_GAMBLR:
    params:
        branch = ", ref = \"" + CFG['inputs']['gamblr_branch'] + "\"" if CFG['inputs']['gamblr_branch'] != "" else "", 
        config_url = CFG["inputs"]["gamblr_config_url"]
    output:
        installed = directory(GAMBLR + "/lib/R/library/GAMBLR"),
        config = "gamblr.yaml"
    conda:
        CFG['conda_envs']['gamblr']
    shell:
        op.as_one_line("""
        wget -qO {output.config} {params.config_url} &&
        R -q -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR"{params.branch})' 
        """)

rule _vcf2maf_deblacklist_maf: 
    input: 
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.raw.maf", 
        gamblr = ancient(rules._vcf2maf_install_GAMBLR.output.installed)
    output: 
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.deblacklisted.maf"
    log: 
        CFG["logs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.deblacklist.log"
    params: 
        blacklist_build = lambda w: VCF2MAF_GENOME_VERSION[w.genome_build], 
        blacklist_template = lambda w: config["lcr-modules"]["vcf2maf"]["inputs"]["blacklist_template"], 
        project_base = CFG["inputs"]["project_base"], 
        seq_type_blacklist = CFG["options"]["seq_type_blacklist"], 
        threshold = CFG["options"]["drop_threshold"] 
    conda:
        CFG['conda_envs']['gamblr']
    threads:
        CFG["threads"]["deblacklist"]
    resources:
        **CFG["resources"]["deblacklist"]
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"]
    script:
        config["lcr-modules"]["vcf2maf"]['inputs']['deblacklisting']
        


def get_add_mafs(wildcards): 
    CFG = config["lcr-modules"]["vcf2maf"]
    # Retreive the patient ID for this tumour
    patient_id = op.filter_samples(CFG["runs"], tumour_sample_id = wildcards.tumour_id).tumour_patient_id.tolist()
    # Retreive all tumours matching this patient
    this_patient = op.filter_samples(CFG["runs"], tumour_patient_id = patient_id, tumour_genome_build = wildcards.genome_build)
    this_patient = this_patient[~(this_patient["tumour_sample_id"].isin([wildcards.tumour_id]) & this_patient["tumour_seq_type"].isin([wildcards.seq_type]))]
    # Subset to one seq type based on config value
    if not CFG["options"]["across_seq_types"]: 
        this_patient = op.filter_samples(this_patient, tumour_seq_type = wildcards.seq_type)
    add_mafs = expand(
        CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.maf", 
        zip, 
        genome_build = this_patient["tumour_genome_build"], 
        seq_type = this_patient["tumour_seq_type"], 
        tumour_id = this_patient["tumour_sample_id"], 
        normal_id = this_patient["normal_sample_id"], 
        pair_status = this_patient["pair_status"], 
        base_name = [CFG["vcf_base_name"]] * len(this_patient.tumour_sample_id),
        allow_missing = True
    )
    return {"add_maf_files": add_mafs}

# Identify the subset of tumours with multiple tumour samples for augmenting 

if CFG["options"]["augment"]: 
    if not CFG["options"]["across_seq_types"]: 
        MULTI_SAMPLES = CFG["runs"]\
            .groupby(["tumour_patient_id", "tumour_seq_type", "tumour_genome_build"])\
            .filter(lambda x: len(x) > 1)\
            .reset_index()
    else: 
        MULTI_SAMPLES = CFG["runs"]\
            .groupby(["tumour_patient_id", "tumour_genome_build"])\
            .filter(lambda x: len(x) > 1)\
            .reset_index()
else: 
    MULTI_SAMPLES = op.discard_samples(CFG["runs"], tumour_sample_id = CFG["runs"]["tumour_sample_id"].tolist())

rule _vcf2maf_augment_maf: 
    input: 
        unpack(get_add_mafs),
        tumour_bam = str(rules._vcf2maf_input_bam.output.bam), 
        index_maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.maf"
    output:
        augmented_maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.augmented_maf"
    params: 
        script = CFG["inputs"]["augment_ssm"]
    threads: 
        CFG['threads']['augment']
    conda: 
        CFG['conda_envs']['pysam']
    resources:
        **CFG['resources']['augment']
    group: "bam_and_augment"
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"], 
        filter = "|".join(["raw", "deblacklisted"]), 
        tumour_id = "|".join(MULTI_SAMPLES.tumour_sample_id.tolist())
    script:
        AUGMENT_SSM


rule _vcf2maf_output_original: 
    input: 
        maf = CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.{maf}"
    output: 
        maf = CFG["dirs"]["outputs"] + "{filter}/{maf}/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{base_name}.maf"
    wildcard_constraints: 
        base_name = CFG["vcf_base_name"], 
        filter = "|".join(["raw", "deblacklisted"]), 
        maf = "|".join(["maf", "augmented_maf"]) 
    run: 
        op.relative_symlink(input.maf, output.maf, in_module = True)


    

def get_original_genome(wildcards):
    # Determine the original (i.e. input) reference genome for this sample
    # Since this module projects to various output genome builds, we need to parse the sample table for the starting build
    # To determine what we need to do
    sample_table = config['lcr-modules']["vcf2maf"]["samples"]
    sample_entry = sample_table.loc[(sample_table["sample_id"] == wildcards.tumour_id) & (sample_table["seq_type"] == wildcards.seq_type)]
    if len(sample_entry) == 0:
        raise AttributeError("Unable to locate a a sample with tumour_id:{wildcards.tumour_id}, normal_id:{wildcards.normal_id}, seq_type:{wildcards.seq_type} in the \'sample\' table")
    original_genome_build = sample_entry.genome_build.tolist()
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
    CFG = config["lcr-modules"]["vcf2maf"]
    original_genome_build = get_original_genome(wildcards)[0]
    return {
        "maf": 
            expand(
                CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.{maf}", 
                **wildcards, 
                genome_build = original_genome_build
            ),
        "convert_coord": config['lcr-modules']["vcf2maf"]["inputs"]["convert_coord"],
        "chains": get_chain(original_genome_build)
    }

rule _vcf2maf_crossmap:
    input:
        unpack(crossmap_input)
    output:
        maf = temp(CFG["dirs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.{filter}.{maf}")
    log:
        stdout = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.{filter}.{maf}.crossmap.stdout.log",
        stderr = CFG["logs"]["crossmap"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.{filter}.{maf}.crossmap.stderr.log"
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


def get_normalize_input(wildcards, todo_only = False):
    CFG = config["lcr-modules"]["vcf2maf"]
    new_genome_build  = wildcards.target_build
    # Since snakemake only knows what the TARGET genome build is, we need to find the source
    original_genome_builds = get_original_genome(wildcards)
    original_genome_version = [VCF2MAF_GENOME_VERSION[x] for x in original_genome_builds]
    new_genome_version = VCF2MAF_GENOME_VERSION[new_genome_build]

    # Do we need to run CrossMap on this? Check the genome version
    if new_genome_build in original_genome_builds:
        # Source matches. CrossMap not necessary. No maf normalization needed. 
        maf = expand(
            CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.{maf}", 
            **wildcards, 
            genome_build = new_genome_build
        )
        todo = "symlink"
    elif new_genome_version in original_genome_version: 
        # Source matches. CrossMap not necessary. 
        # Get the list index of the matching item
        index = original_genome_version.index(new_genome_version)
        # Get the input maf file
        maf = expand(
            CFG["dirs"]["vcf2maf"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.gnomad_filtered.{filter}.{maf}", 
            **wildcards, 
            genome_build = original_genome_builds[index]
        )
        # Check if the chr prefix status is the same e.g. for grch37 vs. hs37d5. 
        # If yes, symlink. If no, normalize. 
        if VCF2MAF_GENOME_PREFIX[new_genome_build] == VCF2MAF_GENOME_PREFIX[original_genome_builds[index]]: 
            todo = "symlink"
        else:
            todo = "normalize"
    else:
        # Source doesn't match. CrossMap and normalization necessary.
        target_build = "hg38" if new_genome_version == "grch38" else "hg19" 
        maf = expand(
            rules._vcf2maf_crossmap.output.maf, 
            target_build = target_build, 
            allow_missing = True
        )
        todo = "normalize"

    if todo_only: 
        return todo
    else: 
        return {"maf": maf}
   

def normalize_prefix(wildcards, input_maf, output_maf, dest_chr):
    maf_open = pd.read_csv(input_maf, sep = "\t")
    # To handle CrossMap weirdness, remove all chr-prefixes and add them back later
    maf_open["Chromosome"] = maf_open["Chromosome"].astype(str).str.replace('chr', '')
    if dest_chr:  # Will evaluate to True if the destination genome is chr-prefixed
        # Add chr prefix
        maf_open['Chromosome'] = 'chr' + maf_open['Chromosome'].astype(str)
    maf_open.to_csv(output_maf, sep="\t", index=False)

# Add or remove chr prefix as necessary
rule _vcf2maf_normalize_prefix:
    input:
        unpack(lambda w: get_normalize_input(w, todo_only = False))
    output:
        maf = CFG["dirs"]["normalize"] + "{seq_type}--{target_build}/{tumour_id}--{normal_id}--{pair_status}/{base_name}.{filter}.{maf}"
    params:
        todo = lambda w: get_normalize_input(w, todo_only = True),
        dest_chr = lambda w: VCF2MAF_GENOME_PREFIX[w.target_build]
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"]), 
        base_name = CFG["vcf_base_name"]
    run:
        if params.todo == "symlink": 
            op.relative_symlink(input.maf, output.maf)
        else: 
            maf_open = pd.read_csv(input.maf[0], sep = "\t")
            print("Normalizing " + input.maf[0])
            # To handle CrossMap weirdness, remove all chr-prefixes and add them back later
            maf_open["Chromosome"] = maf_open["Chromosome"].astype(str).str.replace('chr', '')
            if params.dest_chr:  # Will evaluate to True if the destination genome is chr-prefixed
                # Add chr prefix
                maf_open['Chromosome'] = 'chr' + maf_open['Chromosome'].astype(str)
            print("Writing to " + output.maf)
            maf_open.to_csv(output.maf, sep="\t", index=False)
        

rule _vcf2maf_output_maf:
    input:
        maf = str(rules._vcf2maf_normalize_prefix.output.maf),
    output:
        maf = CFG["dirs"]["outputs"] + "{filter}/{maf}/{seq_type}--projection/{tumour_id}--{normal_id}--{pair_status}.{base_name}.{target_build}.maf"
    wildcard_constraints:
        target_build = "|".join(CFG["options"]["target_builds"]), 
        base_name = CFG["vcf_base_name"]
    run:
        op.relative_symlink(input.maf, output.maf, in_module = True)

print(CFG["options"]["filter"])

# Generates the target sentinels for each run, which generate the symlinks
rule _vcf2maf_all:
    input:
        expand(
            expand(
                str(rules._vcf2maf_output_original.output.maf), 
                zip, 
                seq_type = CFG["runs"]["tumour_seq_type"], 
                genome_build = CFG["runs"]["tumour_genome_build"], 
                tumour_id = CFG["runs"]["tumour_sample_id"], 
                normal_id = CFG["runs"]["normal_sample_id"], 
                pair_status = CFG["runs"]["pair_status"], 
                allow_missing = True
            ),
            maf = "maf", 
            base_name = CFG["vcf_base_name"], 
            filter = CFG["options"]["filter"] 
        ), 
        expand(
            expand(
                str(rules._vcf2maf_output_original.output.maf), 
                zip, 
                seq_type = MULTI_SAMPLES["tumour_seq_type"], 
                genome_build = MULTI_SAMPLES["tumour_genome_build"], 
                tumour_id = MULTI_SAMPLES["tumour_sample_id"], 
                normal_id = MULTI_SAMPLES["normal_sample_id"], 
                pair_status = MULTI_SAMPLES["pair_status"], 
                allow_missing = True
            ),
            base_name = CFG["vcf_base_name"], 
            maf = "augmented_maf", 
            filter = CFG["options"]["filter"] 
        ), 
        expand(
            expand(
                str(rules._vcf2maf_output_maf.output.maf), 
                zip,
                seq_type = CFG["runs"]["tumour_seq_type"],
                tumour_id = CFG["runs"]["tumour_sample_id"],
                normal_id = CFG["runs"]["normal_sample_id"],
                pair_status = CFG["runs"]["pair_status"],
                base_name = [CFG["vcf_base_name"]] * len(CFG["runs"]["tumour_sample_id"]),
                allow_missing = True
            ),
            target_build = CFG["options"]["target_builds"], 
            base_name = CFG["vcf_base_name"], 
            maf = "maf", 
            filter = CFG["options"]["filter"] 
        ), 
        expand(
            expand(
                str(rules._vcf2maf_output_maf.output.maf), 
                zip,
                seq_type = MULTI_SAMPLES["tumour_seq_type"],
                tumour_id = MULTI_SAMPLES["tumour_sample_id"],
                normal_id = MULTI_SAMPLES["normal_sample_id"],
                pair_status = MULTI_SAMPLES["pair_status"],
                base_name = [CFG["vcf_base_name"]] * len(MULTI_SAMPLES["tumour_sample_id"]),
                allow_missing = True
            ),
            target_build = CFG["options"]["target_builds"], 
            base_name = CFG["vcf_base_name"], 
            maf = "augmented_maf", 
            filter = CFG["options"]["filter"]
        )
        # Why are there two expand statements? Well we want every iteration of these MAFs for all the target genome builds
        # But that is the only wildcard we want to expand to all iterations
        # Hence the inner expand is using zip()

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
