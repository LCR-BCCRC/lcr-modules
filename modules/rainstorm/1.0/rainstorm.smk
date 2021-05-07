#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Kostiantyn Dreval
# Module Author:    Kostiantyn Dreval
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
# `CFG` is a shortcut to `config["lcr-modules"]["rainstorm"]`
CFG = op.setup_module(
    name = "rainstorm",
    version = "1.0",
    subdirectories = ["inputs", "rainstorm", "doppler", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _rainstorm_install,
    _rainstorm_input_maf,
    _rainstorm_subset_index,
    _rainstorm_output_tsv,
    _rainstorm_all


##### RULES #####

# download rainstorm files from git repo without cloning the repo itself
# decompress files into the 00-inputs
rule _rainstorm_install:
    output:
        rainstorm = CFG["dirs"]["inputs"] + "mutation_rainstorm-0.3/rainstorm.py", # one of the files in the repo
        peaks = CFG["dirs"]["inputs"] + "mutation_rainstorm-0.3/rainstorm_peaks.R" # another script from the repo
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
        maf = CFG["inputs"]["sample_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{cohort_name}.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)


# subset index to only standard chromosomes, so rainstorm will not look for mutations in the small non-standard contigs
# currently, this excludes chromosome Y
rule _rainstorm_subset_index:
    input:
        index = reference_files("genomes/{genome_build}/genome_fasta/genome.fa.fai")
    output:
        index_subset = CFG["dirs"]["inputs"] + "index_files/index_subset.{genome_build}.fa.fai"
    wildcard_constraints:
        genome_build = CFG["options"]["genome_build"]
    shell:
        op.as_one_line("""
        cat {input.index} | perl -ne 'print if /^#|^(chr)*[\dX]+\s.+/' >
        {output.index_subset}
        """)

# main rule to run rainstorm analysis
rule _rainstorm_run:
    input:
        rainstorm = str(rules._rainstorm_install.output.rainstorm),    
        maf = str(rules._rainstorm_input_maf.output.maf),
        index_subset = str(rules._rainstorm_subset_index.output.index_subset)
    output:
        complete = CFG["dirs"]["rainstorm"] + "{genome_build}/{cohort_name}_out_background_100k_binned_density.tsv"
    log:
        stdout = CFG["logs"]["rainstorm"] + "{genome_build}/{cohort_name}_rainstorm.stdout.log",
        stderr = CFG["logs"]["rainstorm"] + "{genome_build}/{cohort_name}_rainstorm.stderr.log"
    params:
        out_name = CFG["dirs"]["rainstorm"] + "{genome_build}/{cohort_name}_out",
        rainstorm_flags = CFG["options"]["flags"]
    conda:
        CFG["conda_envs"]["rainstorm"]
    threads:
        CFG["threads"]["rainstorm"]
    resources:
        **CFG["resources"]["rainstorm"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.cohort_name} on $(hostname) at $(date)" >> {log.stdout};
        python3 {input.rainstorm} maf {input.maf}
        --output_base_name {params.out_name}
        --cpu_num {threads}
        --genome_fai {input.index_subset}
        {params.rainstorm_flags}
        >> {log.stdout}
        2>> {log.stderr} &&
        touch {output.complete}
        echo "DONE {rule} for {wildcards.cohort_name} on $(hostname) at $(date)" >> {log.stdout};
        """)


# generate the list of chromosomes to expand for doppler run
chroms = list(map(str, range(1, 22))) + ["X"] 
chroms = ["chr" + chrom for chrom in chroms]


# Run doppler for the population-based analysis
rule _rainstorm_run_doppler:
    input:
        peaks = str(rules._rainstorm_install.output.peaks),
        tsv = str(rules._rainstorm_run.output.complete),
        maf = str(rules._rainstorm_input_maf.output.maf)
    output:
        tsv = CFG["dirs"]["doppler"] + "{genome_build}/{cohort_name}_mean_waveletSummary_withMaf.tsv"
    log:
        stdout = CFG["logs"]["doppler"] + "{genome_build}/{cohort_name}_rainstorm.stdout.log",
        stderr = CFG["logs"]["doppler"] + "{genome_build}/{cohort_name}_rainstorm.stderr.log"
    params:
        tsv = expand(str(CFG["dirs"]["rainstorm"] + "{{genome_build}}/{{cohort_name}}_out_rainstorm_k_4_mean_{chromosome}.tsv"), zip, chromosome=chroms),
        str_split = CFG["dirs"]["rainstorm"] + "{genome_build}/{cohort_name}_out_rainstorm_k_4_mean_",
        out_name = CFG["dirs"]["doppler"] + "{genome_build}/" + CFG["options"]["cohort_name"] + "_mean_"
    conda:
        CFG["conda_envs"]["rainstorm"]
    threads:
        CFG["threads"]["rainstorm"]
    resources:
        **CFG["resources"]["rainstorm"]
    shell:
        op.as_one_line("""
        echo "running {rule} for {wildcards.cohort_name} on $(hostname) at $(date)" >> {log.stdout};
        Rscript {input.peaks} {params.tsv}
        --stringSplit {params.str_split}
        --input_maf {input.maf}
        --output_base_file {params.out_name}
        >> {log.stdout}
        2>> {log.stderr} &&
        echo "DONE {rule} for {wildcards.cohort_name} on $(hostname) at $(date)" >> {log.stdout};
        """)


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _rainstorm_output_tsv:
    input:
        tsv = str(rules._rainstorm_run_doppler.output.tsv)
    output:
        tsv = CFG["dirs"]["outputs"] + "tsv/{genome_build}/{cohort_name}.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _rainstorm_all:
    input:
        expand(
            [
                str(rules._rainstorm_output_tsv.output.tsv),
            ],
            zip,  # Run expand() with zip(), not product()
            genome_build=CFG["options"]["genome_build"],
            cohort_name=CFG["options"]["cohort_name"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
