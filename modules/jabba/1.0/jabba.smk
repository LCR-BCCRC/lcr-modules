#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
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
# `CFG` is a shortcut to `config["lcr-modules"]["jabba"]`
CFG = op.setup_module(
    name = "jabba",
    version = "1.0",
    subdirectories = ["inputs", "fragcounter", "dryclean", "outputs"],
)


# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _jabba_install_fragcounter,
    _jabba_install_dryclean,
    _jabba_input_bam,
    _jabba_link_normal_rds,
    _jabba_link_dryclean_normal_rds,
    _jabba_link_dryclean_tumour_rds,
    _jabba_all


##### RULES #####

### Adjust coverage for GC and mappability ###

# Install JaBba suite from github
rule _jabba_install_fragcounter:
    output:
        complete = CFG["dirs"]["fragcounter"] + "fragcounter.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'remotes::install_github("mskilab/fragCounter")' 
            &&
        touch {output.complete}
        """)

rule _jabba_install_dryclean:
    output:
        complete = CFG["dirs"]["dryclean"] + "dryclean.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)' 
                -e 'remotes::install_github("mskilab/dryclean")'
            &&
        touch {output.complete}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _jabba_input_bam:
    input:
        bam = CFG["inputs"]["sample_bam"],
        bai = CFG["inputs"]["sample_bai"]
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{sample_id}.bai",
        crai = CFG["dirs"]["inputs"]+ "bam/{seq_type}--{genome_build}/{sample_id}.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)

# Runs fragcounter on individual samples
rule _jabba_run_fragcounter:
    input:
        installed = str(rules._jabba_install_fragcounter.output.complete),
        bam = str(rules._jabba_input_bam.output.bam),
        gc = reference_files("genomes/{genome_build}/annotations/jabba/gc1000.rds"),
        map = reference_files("genomes/{genome_build}/annotations/jabba/map1000.rds")
    output:
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/cov.rds"
    log:
        stdout = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/fc.stdout.log",
        stderr = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{sample_id}/fc.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["fragcounter"]
    resources:
        mem_mb = CFG["mem_mb"]["fragcounter"]
    shell:
        op.as_one_line(""" 
        export PATH=$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))'):$PATH
            &&
        frag -b {input.bam} -d `dirname {input.gc}` -w 1000 -o `dirname {output.rds}` > {log.stdout} 2> {log.stderr}
        """)


### Make panel-of-normal and germline filter ###


# Symlink fragcounter outputs of normals; required for PON and germline.filter file generation
rule _jabba_link_normal_rds:
    input:
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{normal_id}/cov.rds"
    output:
        rds = CFG["dirs"]["fragcounter"] + "cov/{seq_type}--{genome_build}/normal/{normal_id}.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Generate Panel-of-normals from fragcounter normal outputs
rule _jabba_make_pon:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        par = reference_files("genomes/{genome_build}/annotations/jabba/PAR_{genome_build}.rds"),
        rds = expand(CFG["dirs"]["fragcounter"] + "cov/{{seq_type}}--{{genome_build}}/normal/{normal_id}.cov.rds", normal_id = CFG["runs"]["normal_sample_id"]),
        make_pon = CFG["inputs"]["make_pon"]
    output:
        tbl = CFG["dirs"]["dryclean"] + "pon/{seq_type}--{genome_build}/normal_table.rds",
        pon = CFG["dirs"]["dryclean"] + "pon/{seq_type}--{genome_build}/detergent.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "pon/{seq_type}--{genome_build}/make_pon.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "pon/{seq_type}--{genome_build}/make_pon.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["make_pon"]
    resources:
        mem_mb = CFG["mem_mb"]["make_pon"]
    shell:
        op.as_one_line("""
        Rscript {input.make_pon} {threads} `dirname {input.rds[0]}` {output.tbl} {output.pon} {input.par} {wildcards.genome_build} > {log.stdout} 2> {log.stderr}
        """)

# Run dryclean decomposition on normals; used to identify germline events and artifacts
rule _jabba_run_dryclean_normal:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        rds = str(rules._jabba_run_fragcounter.output.rds),
        pon = str(rules._jabba_make_pon.output.pon)
    output:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}/drycleaned.cov.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["dryclean"]
    resources:
       mem_mb = CFG["mem_mb"]["dryclean"]
    wildcard_constraints:
        sample_id = "|".join(CFG["runs"]["normal_sample_id"])
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, is.human = FALSE, mc.cores = {threads})'
                -e 'saveRDS(decomp, "{output.rds}")' > {log.stdout} 2> {log.stderr}
        """)

# Symlink dryclean normal outputs; needed for generating germline.filter file
rule _jabba_link_dryclean_normal_rds:
    input:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{normal_id}/drycleaned.cov.rds"
    output:
        rds = CFG["dirs"]["dryclean"] + "cov/{seq_type}--{genome_build}/normal/{normal_id}.drycleaned.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Make germline filter file from drycleaned normal outputs
rule _jabba_make_germline_filter:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        tbl = str(rules._jabba_make_pon.output.tbl),
        rds = expand(CFG["dirs"]["dryclean"] + "cov/{{seq_type}}--{{genome_build}}/normal/{normal_id}.drycleaned.cov.rds", normal_id = CFG["runs"]["normal_sample_id"]),
        make_germline = CFG["inputs"]["make_germline"]
    output:
        germline = CFG["dirs"]["dryclean"] + "germline_filter/{seq_type}--{genome_build}/germline.markers.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "germline_filter/{seq_type}--{genome_build}/germline_filter.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "germline_filter/{seq_type}--{genome_build}/germline_filter.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["make_germline"]
    resources:
        mem_mb = CFG["mem_mb"]["make_germline"]
    shell:
        op.as_one_line(""" 
        Rscript {input.make_germline} {input.tbl} `dirname {input.rds[0]}` {output.germline} 0.5 0.98 > {log.stdout} 2> {log.stderr}
        """)


### Run dryclean on tumours to remove background and germline variation ###


rule _jabba_run_dryclean_tumour:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        rds = str(rules._jabba_run_fragcounter.output.rds),
        pon = str(rules._jabba_make_pon.output.pon),
        germline = str(rules._jabba_make_germline_filter.output.germline)
    output:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}/drycleaned.cov.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{sample_id}.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["dryclean"]
    resources:
        mem_mb = CFG["mem_mb"]["dryclean"]
    wildcard_constraints:
        sample_id = "|".join(CFG["runs"]["tumour_sample_id"])
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, is.human = FALSE, mc.cores = {threads}, germline.file = "{input.germline}")'
                -e 'saveRDS(decomp, "{output.rds}")' > {log.stdout} 2> {log.stderr}
        """)

rule _jabba_link_dryclean_tumour_rds:
    input:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}/drycleaned.cov.rds"
    output:
        rds = CFG["dirs"]["outputs"] + "rds/{seq_type}--{genome_build}/{tumour_id}.drycleaned.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Generates the target sentinels for each run, which generate the symlinks
rule _jabba_all:
    input:
        expand(
            [
                str(rules._jabba_link_dryclean_tumour_rds.output.rds),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
