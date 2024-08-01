#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
# Contributors:     Sierra Gillis


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
from datetime import datetime
import numpy as np

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
# `CFG` is a shortcut to `config["lcr-modules"]["rainstorm"]`
CFG = op.setup_module(
    name = "rainstorm",
    version = "1.1",
    subdirectories = ["inputs", "rainstorm", "doppler", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _rainstorm_install,
    _rainstorm_input_maf,
    _rainstorm_input_subsetting_categories,
    _rainstorm_subset_index,
    _rainstorm_prepare_maf,
    _rainstorm2bed,
    _rainstorm_output,
    _rainstorm_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS =  os.path.abspath(config["lcr-modules"]["rainstorm"]["prepare_mafs"])

# Download rainstorm files from git repo without cloning the repo itself, decompress files into the 00-inputs
rule _rainstorm_install:
    output:
        rainstorm = CFG["dirs"]["inputs"] + "mutation_rainstorm-" + str(CFG["options"]["rainstorm_version"]) + "/rainstorm.py", # one of the files in the repo
        peaks = CFG["dirs"]["inputs"] + "mutation_rainstorm-" + str(CFG["options"]["rainstorm_version"]) + "/rainstorm_peaks.R", # another script from the repo
        doppler2bed = CFG["dirs"]["inputs"] + "mutation_rainstorm-" + str(CFG["options"]["rainstorm_version"]) + "/doppler_to_bed.py" # convert output to bed format
    params:
        url = "https://github.com/rdmorin/mutation_rainstorm/archive/refs/tags/v" + str(CFG["options"]["rainstorm_version"]) + ".tar.gz",
        folder = CFG["dirs"]["inputs"]
    shell:
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _rainstorm_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{sample_set}--{genome_build}--{launch_date}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

# Subset index to only standard chromosomes, so rainstorm will not look for mutations in the small non-standard contigs
# currently, this excludes chromosome Y
rule _rainstorm_subset_index:
    input:
        index = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        index_subset = CFG["dirs"]["inputs"] + "index_files/index_subset.{genome_build}.fa.fai"
    shell:
        op.as_one_line("""
        cat {input.index} | perl -ne 'print if /^#|^(chr)*[\dX]+\s.+/' >
        {output.index_subset}
        """)

# Symlinks the subsetting categories file into the module results directory (under '00-inputs/')
rule _rainstorm_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

# Prepare the maf file for the input to rainstorm
checkpoint _rainstorm_prepare_maf:
    input:
        maf = str(rules._rainstorm_input_maf.output.maf),
        subsetting_categories = str(rules._rainstorm_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/done"
    log:
        CFG["logs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "rainstorm",
        metadata_cols = CFG["samples"],
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

# Main rule to run rainstorm analysis
rule _rainstorm_run:
    input:
        rainstorm = str(rules._rainstorm_install.output.rainstorm),
        maf = CFG["dirs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/{md5sum}.maf",
        content = CFG["dirs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/{md5sum}.maf.content",
        index_subset = str(rules._rainstorm_subset_index.output.index_subset)
    output:
        complete = CFG["dirs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_out_background_100k_binned_density.tsv"
    log:
        stdout = CFG["logs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_rainstorm.stdout.log",
        stderr = CFG["logs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_rainstorm.stderr.log"
    params:
        out_name = CFG["dirs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_out",
        rainstorm_flags = CFG["options"]["flags"]
    conda:
        CFG["conda_envs"]["rainstorm"]
    threads:
        CFG["threads"]["rainstorm"]
    resources:
        **CFG["resources"]["rainstorm"]
    shell:
        op.as_one_line("""
        python3 {input.rainstorm} maf {input.maf}
        --output_base_name {params.out_name}
        --cpu_num {threads}
        --genome_fai {input.index_subset}
        {params.rainstorm_flags}
        >> {log.stdout}
        2>> {log.stderr} &&
        touch {output.complete}
        """)

# Generate the list of chromosome files for doppler run
def _get_chrom_tsvs(wildcards):
    CFG = config["lcr-modules"]["rainstorm"]

    if wildcards.genome_build == "hg38":
        chroms = list(map(str, range(1, 22))) + ["X"]
        chroms = ["chr" + chrom for chrom in chroms]
    elif wildcards.genome_build == "grch37":
        chroms = list(map(str, range(1, 22))) + ["X"]

    return expand(
        CFG["dirs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_out_rainstorm_k_4_mean_{chromosome}.tsv",
        genome_build = wildcards.genome_build,
        sample_set = wildcards.sample_set,
        launch_date = wildcards.launch_date,
        md5sum = wildcards.md5sum,
        chromosome=chroms
        )

# Run doppler for the population-based analysis
rule _rainstorm_run_doppler:
    input:
        peaks = str(rules._rainstorm_install.output.peaks),
        tsv = str(rules._rainstorm_run.output.complete),
        maf = CFG["dirs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/{md5sum}.maf"
    output:
        tsv = CFG["dirs"]["doppler"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_mean_waveletSummary_withMaf.tsv"
    log:
        stdout = CFG["logs"]["doppler"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_rainstorm.stdout.log",
        stderr = CFG["logs"]["doppler"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_rainstorm.stderr.log"
    params:
        tsv = _get_chrom_tsvs,
        str_split = CFG["dirs"]["rainstorm"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_out_rainstorm_k_4_mean_",
        out_name = CFG["dirs"]["doppler"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}" + "_mean_"
    conda:
        CFG["conda_envs"]["rainstorm"]
    threads:
        CFG["threads"]["rainstorm"]
    resources:
        **CFG["resources"]["rainstorm"]
    shell:
        op.as_one_line("""
        Rscript --vanilla {input.peaks} {params.tsv}
        --stringSplit {params.str_split}
        --input_maf {input.maf}
        --output_base_file {params.out_name}
        >> {log.stdout}
        2>> {log.stderr}
        """)

# Convert doppler output to bed file
rule _rainstorm2bed:
    input:
        converter = str(rules._rainstorm_install.output.doppler2bed),
        doppler = str(rules._rainstorm_run_doppler.output.tsv)
    output:
        bed = CFG["dirs"]["doppler"] + "{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}_mean_waveletSummary_withMaf.bed"
    conda:
        CFG["conda_envs"]["rainstorm"]
    shell:
        op.as_one_line("""
        python3 {input.converter}
        {input.doppler}
        -o {output.bed}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _rainstorm_output:
    input:
        tsv = str(rules._rainstorm_run_doppler.output.tsv),
        bed = str(rules._rainstorm2bed.output.bed)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}.tsv",
        bed = CFG["dirs"]["outputs"] + "bed/{genome_build}/{sample_set}--{launch_date}--{md5sum}/{sample_set}.bed"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)
        op.relative_symlink(input.bed, output.bed, in_module= True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["rainstorm"]
    checkpoint_output = os.path.dirname(str(checkpoints._rainstorm_prepare_maf.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum}.maf.content")
    return expand(
        [
            CFG["dirs"]["outputs"] + "tsv/{{genome_build}}/{{sample_set}}--{{launch_date}}--{md5sum}/{{sample_set}}.tsv",
            CFG["dirs"]["outputs"] + "bed/{{genome_build}}/{{sample_set}}--{{launch_date}}--{md5sum}/{{sample_set}}.bed"
        ],
        md5sum = SUMS
        )

# Aggregates outputs to remove md5sum from rule all
rule _rainstorm_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{sample_set}--{genome_build}--{launch_date}.done"
    shell:
        op.as_one_line("""touch {output.aggregate}""")

# Generates the target sentinels for each run, which generate the symlinks
rule _rainstorm_all:
    input:
        expand(
            [
                CFG["dirs"]["inputs"] + "{sample_set}--{genome_build}--{launch_date}/done",
                str(rules._rainstorm_aggregate.output.aggregate)
            ], # no "zip" so that genome builds are expanded correctly if multiple are given
            genome_build=CFG["options"]["genome_build"],
            sample_set=CFG["sample_set"],
            launch_date = launch_date)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
