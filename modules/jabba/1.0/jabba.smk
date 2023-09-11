#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     Krysta Coyle


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

##make normals table and get rid of tissue_status

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["jabba"]`
CFG = op.setup_module(
    name = "jabba",
    version = "1.0",
    subdirectories = ["inputs","setup","pon","fragcounter", "dryclean", "jabba",  "outputs"],
)

os.environ["CPLEX_DIR"] = str(CFG["CPLEX_DIR"])

# Define rules to be run locally when using a compute cluster
localrules:
    _jabba_install_dryclean,
    _jabba_install_fragcounter,
    _jabba_pon_symlink_normal_bams,
    _jabba_pon_link_dryclean_normal_rds,
    _jabba_pon_symlink_fragcounter,
    _jabba_install_jabba,
    _jabba_input_tumour_bam,
    _jabba_link_graph_rds,
    _jabba_all

##### RULES #####

#Creating table of normals because wildcard_constraints did not work as expected
NORMALS = op.filter_samples(SAMPLES, tissue_status="normal")

# Install JaBba suite from github
rule _jabba_install_fragcounter:
    output:
        complete = CFG["dirs"]["setup"] + "fragcounter.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'Sys.unsetenv("GITHUB_PAT")'
                -e 'if ("fragCounter" %in% installed.packages()==FALSE) {{remotes::install_github("morinlab/fragCounter", upgrade = TRUE, force = TRUE)}}'
                -e 'library(fragCounter)'
            &&
        touch {output.complete}
        """)

rule _jabba_install_dryclean:
    input:
        installed1 = str(rules._jabba_install_fragcounter.output.complete)
    output:
        complete = CFG["dirs"]["setup"] + "dryclean.installed"
    conda: CFG["conda_envs"]["jabba"]
    shell:
        op.as_one_line("""
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'Sys.unsetenv("GITHUB_PAT")'
                -e 'if ("dryclean" %in% installed.packages()==FALSE) {{remotes::install_github("morinlab/dryclean", upgrade = TRUE, force = TRUE)}}'
                -e 'library(dryclean)'
            &&
        touch {output.complete}
        """)


rule _jabba_install_jabba:
    input:
        installed1 = str(rules._jabba_install_dryclean.output.complete),
        installed2 = str(rules._jabba_install_fragcounter.output.complete)
    output:
        complete = CFG["dirs"]["jabba"] + "jabba.installed"
    conda: CFG["conda_envs"]["jabba"]
    params:
        cplex_dir = CFG["CPLEX_DIR"]
    shell:
        op.as_one_line(""" 
        Rscript -e 'Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)'
                -e 'Sys.setenv(CPLEX_DIR = "{params.cplex_dir}")'
                -e 'Sys.unsetenv("GITHUB_PAT")'
                -e 'if ("copynumber" %in% installed.packages()==FALSE) {{BiocManager::install("copynumber")}}'
                -e 'if (is.null(packageDescription("copynumber")$GithubUsername)) {{remotes::install_github("ShixiangWang/copynumber", dependencies = TRUE)}}'
                -e 'if ("JaBbA" %in% installed.packages()==FALSE) {{remotes::install_github("mskilab/JaBbA", upgrade = TRUE, force = TRUE)}}'
                -e 'library(JaBbA)'
            &&
        touch {output.complete}
        """)


##### JaBbA Panel-of-Normals (PON) generation ######

rule _jabba_pon_symlink_normal_bams:
    input: 
        bam = "data/{seq_type}_bams/{sample_id}.{genome_build}.bam",
        bai = "data/{seq_type}_bams/{sample_id}.{genome_build}.bam.bai"
    output: 
        bam = CFG["dirs"]["pon"] + "00-normal_bams/{seq_type}--{genome_build}/{tissue_status}/{sample_id}.bam",
        bai = CFG["dirs"]["pon"] + "00-normal_bams/{seq_type}--{genome_build}/{tissue_status}/{sample_id}.bam.bai",
        crai = CFG["dirs"]["pon"] + "00-normal_bams/{seq_type}--{genome_build}/{tissue_status}/{sample_id}.bam.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)

# Run fragcounter on normals to get GC/mappability corrected coverage values
rule _jabba_pon_run_fragcounter:
    input:
        installed = str(rules._jabba_install_fragcounter.output.complete),
        bam = CFG["dirs"]["pon"] + "00-normal_bams/{seq_type}--{genome_build}/{tissue_status}/{sample_id}.bam",
        gc = reference_files("genomes/{genome_build}/annotations/jabba/gc1000.rds"),
        map = reference_files("genomes/{genome_build}/annotations/jabba/map1000.rds"),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        rds = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/cov.rds"
    log:
        stdout = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/stdout.txt",
        stderr = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/stderr.txt",
    params:
        dir = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}"
    conda: CFG["conda_envs"]["jabba"]
    threads: 1
    # wildcard_constraints:
    #     tissue_status = "normal"
    resources:
        mem_mb = 4096
    shell:
        op.as_one_line("""
        unset R_HOME &&
        mkdir -p {params.dir} &&
        FRAG=$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))'); Rscript $FRAG/frag -b {input.bam} -r {input.ref} -w 1000 -d `dirname {input.gc}` -o `dirname {output.rds}` > {log.stdout} 2> {log.stderr}
        """)

rule _jabba_pon_symlink_fragcounter:
    input:
        rds = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/cov.rds"
    output:
        rds = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/cov/{sample_id}.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Use dryclean to create panel-of-normal (PON) from fragcounter coverages
rule _jabba_pon_make_pon:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        par = reference_files("genomes/{genome_build}/annotations/jabba/PAR_{genome_build}.rds"),
        make_pon = CFG["inputs"]["make_pon"],
        rds = expand(
            [
                str(rules._jabba_pon_symlink_fragcounter.output.rds)
            ],
            zip,
            genome_build=NORMALS["genome_build"],
            sample_id=NORMALS["sample_id"],
            seq_type=NORMALS["seq_type"],
            tissue_status= NORMALS["tissue_status"])
    output:
        tbl = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{tissue_status}/normal_table.rds",
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{tissue_status}/detergent.rds"
    params:
        choose_samples = 'cluster',
        dir = CFG["dirs"]["pon"] + "01-fragcounter/{seq_type}--{genome_build}/{tissue_status}/cov/"
    conda: CFG["conda_envs"]["jabba"]
    threads: 10
    resources:
       mem_mb = 200000 #may need to update if number of normals is >250
    shell:
        op.as_one_line("""
        Rscript {input.make_pon} {threads} {params.dir} {output.tbl} {output.pon} {input.par} {wildcards.genome_build} {params.choose_samples}
        """)

# Improve coverage signal fidelity by running each normal against PON
# Serves to identify and separate background variations (noise) from foreground variation (signal)
rule _jabba_pon_run_dryclean_normal:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        rds = str(rules._jabba_pon_run_fragcounter.output.rds),
        pon = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{tissue_status}/detergent.rds"
    output:
        rds = CFG["dirs"]["pon"] + "02-dryclean/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/drycleaned.cov.rds"
    conda: CFG["conda_envs"]["jabba"]
    threads: 12
    resources:
       mem_mb = 30000
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, mc.cores = {threads}, germline.filter = FALSE)'
                -e 'saveRDS(decomp, "{output.rds}")'
        """)


rule _jabba_pon_link_dryclean_normal_rds:
    input:
        rds = CFG["dirs"]["pon"] + "02-dryclean/{seq_type}--{genome_build}/{tissue_status}/run/{sample_id}/drycleaned.cov.rds"
    output:
        rds = CFG["dirs"]["pon"] + "02-dryclean/{seq_type}--{genome_build}/{tissue_status}/cov/{sample_id}.drycleaned.cov.rds"
    run:
        op.relative_symlink(input.rds, output.rds)

# Use dryclean normals to identify germline regions which will be used
# on tumour data to find germline regions.
# Save normals to new directory!
rule _jabba_pon_make_germline_filter:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        tbl = str(rules._jabba_pon_make_pon.output.tbl),
        rds = expand(
            [
                str(rules._jabba_pon_link_dryclean_normal_rds.output.rds)
            ],
            zip,
            genome_build=NORMALS["genome_build"],
            sample_id=NORMALS["sample_id"],
            tissue_status=NORMALS["tissue_status"],
            seq_type=NORMALS["seq_type"]
            ),
        make_germline = CFG["inputs"]["make_germline"]
    output:
        germline = CFG["dirs"]["pon"] + "{seq_type}--{genome_build}/{tissue_status}/germline.markers.rds"
    conda: CFG["conda_envs"]["jabba"]
    threads: 25
    resources:
        mem_mb = 300000
    params:
        dir = CFG["dirs"]["pon"] + "02-dryclean/{seq_type}--{genome_build}/{tissue_status}/cov/"
    shell:
        op.as_one_line("""
        Rscript {input.make_germline} {input.tbl} {params.dir} {output.germline} 0.5 0.1 {threads}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _jabba_input_tumour_bam:
    input:
        bam = "data/{seq_type}_bams/{tumour_id}.{genome_build}.bam",
        bai = "data/{seq_type}_bams/{tumour_id}.{genome_build}.bam.bai"
    output:
        bam = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bam",
        bai = CFG["dirs"]["inputs"] + "bam/{seq_type}--{genome_build}/{tumour_id}.bai",
        crai = CFG["dirs"]["inputs"]+ "bam/{seq_type}--{genome_build}/{tumour_id}.crai"
    run:
        op.relative_symlink(input.bam, output.bam)
        op.relative_symlink(input.bai, output.bai)
        op.relative_symlink(input.bai, output.crai)


rule _jabba_input_gridss_vcf:
    input:
        junc = CFG["inputs"]["gridss_junc"]
    output:
        junc = CFG["dirs"]["inputs"] + "junc/gridss/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf"
    run:
        op.relative_symlink(input.junc, output.junc)

rule _jabba_process_gridss_vcf:
    input:
        junc = str(rules._jabba_input_gridss_vcf.output.junc)
    output:
        junc = CFG["dirs"]["inputs"] + "junc/gridss/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bnd.vcf"
    shell:
        op.as_one_line("""
        zcat {input.junc} | awk 'BEGIN {{OFS=FS=\"\t\"}} $0 ~ /^#/ || $3 !~ /b$/' > {output.junc}
        """)

rule _jabba_input_manta_vcf:
    input:
        junc = CFG["inputs"]["manta_junc"]
    output:
        junc = CFG["dirs"]["inputs"] + "junc/manta/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.vcf"
    run:
        op.relative_symlink(input.junc, output.junc)


rule _jabba_merge_svs:
    input:
        installed = str(rules._jabba_install_jabba.output.complete),
        manta = str(rules._jabba_input_manta_vcf.output.junc),
        gridss = str(rules._jabba_process_gridss_vcf.output.junc)
    output:
        junc = CFG["dirs"]["inputs"] + "junc/merged/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.rds"
    log:
        stdout = CFG["logs"]["inputs"] + "run/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merge_svs.stdout.log",
        stderr = CFG["logs"]["inputs"] + "run/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/merge_svs.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["merge_svs"]
    params:
        merge_svs = CFG["inputs"]["merge_svs"]
    resources:
        mem_mb = CFG["mem_mb"]["merge_svs"]
    shell:
        op.as_one_line("""
        Rscript {params.merge_svs} --manta {input.manta} --gridss {input.gridss} --genome {wildcards.genome_build} --rds {output} > {log.stdout} 2> {log.stderr}
        """)


# Runs fragcounter on individual samples
rule _jabba_run_fragcounter:
    input:
        installed = str(rules._jabba_install_fragcounter.output.complete),
        bam = str(rules._jabba_input_tumour_bam.output.bam),
        gc = reference_files("genomes/{genome_build}/annotations/jabba/gc1000.rds"),
        map = reference_files("genomes/{genome_build}/annotations/jabba/map1000.rds"),
        ref = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{tumour_id}/cov.rds"
    log:
        stdout = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{tumour_id}/fc.stdout.log",
        stderr = CFG["logs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{tumour_id}/fc.stderr.log"
    wildcard_constraints:
        pair_status = "matched"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["fragcounter"]
    resources:
        mem_mb = CFG["mem_mb"]["fragcounter"]
    shell:
        op.as_one_line(""" 
        FRAG=$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))'); $FRAG/frag -b {input.bam} -r {input.ref} -w 1000 -d `dirname {input.gc}` -o `dirname {output.rds}` > {log.stdout} 2> {log.stderr}
        """)


### Run dryclean on tumours to remove background and germline variation ###

rule _jabba_run_dryclean_tumour:
    input:
        installed = str(rules._jabba_install_dryclean.output.complete),
        rds = CFG["dirs"]["fragcounter"] + "run/{seq_type}--{genome_build}/{tumour_id}/cov.rds",
        germline = expand(
            [
                str(rules._jabba_pon_make_germline_filter.output.germline)
            ],
            zip,
            tissue_status=str("normal"),
            allow_missing=True),
        pon = expand(
            [
                str(rules._jabba_pon_make_pon.output.pon)
            ],
            zip,
            tissue_status=str("normal"),
            allow_missing=True)
    output:
        rds = CFG["dirs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}/drycleaned.cov.rds"
    log:
        stdout = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}.stdout.log",
        stderr = CFG["logs"]["dryclean"] + "run/{seq_type}--{genome_build}/{tumour_id}.stderr.log"
    conda: CFG["conda_envs"]["jabba"]
    threads: CFG["threads"]["dryclean"]
    resources:
        mem_mb = CFG["mem_mb"]["dryclean"]
    wildcard_constraints:
        tumour_id = "|".join(CFG["runs"]["tumour_sample_id"])
    shell:
        op.as_one_line("""
        Rscript -e 'library(dryclean); library(parallel)'
                -e 'samp <- readRDS("{input.rds}")'
                -e 'decomp <- start_wash_cycle(cov = samp, detergent.pon.path = "{input.pon}", whole_genome = TRUE, mc.cores = {threads}, germline.file = "{input.germline}", germline.filter = TRUE)'
                -e 'decomp <- gUtils::gr2dt(decomp)'
                -e 'saveRDS(decomp, "{output.rds}")' > {log.stdout} 2> {log.stderr}
        """)

def get_pp(wildcards):
    pp = config['lcr-modules']['jabba']['pp']

    purity = pp[wildcards.genome_build][wildcards.tumour_id]['purity']
    purity_range = [round(purity - 0.1, 2), round(purity + 0.1, 2)]
    purity_range = [str(1) if x > 1 else str(0) if x < 0 else str(x) for x in purity_range]

    ploidy = pp[wildcards.genome_build][wildcards.tumour_id]['ploidy']

    return({'purity': ','.join(purity_range), 'ploidy': ploidy})

rule _jabba_run_jabba:
    input:
        installed = str(rules._jabba_install_jabba.output.complete),
        rds = str(rules._jabba_run_dryclean_tumour.output.rds),
        junc = str(rules._jabba_merge_svs.output.junc)
    output:
        rds = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/jabba.simple.gg.rds",
        karyograph = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/karyograph.rds",
        simple = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/jabba.simple.rds"
    log:
        stdout = CFG["logs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stdout.log",
        stderr = CFG["logs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.stderr.log"
    params:
        opts = CFG["options"]["jabba"],
        purity = lambda w: get_pp(w)['purity'],
        ploidy = lambda w: get_pp(w)['ploidy'],#--purity {params.purity} 
        CPLEX = CFG["CPLEX_DIR"]
    conda: CFG["conda_envs"]["jabba"]
    wildcard_constraints:
        pair_status = "matched"
    threads: CFG["threads"]["jabba"]
    resources:
        mem_mb = CFG["mem_mb"]["jabba"],
        license = 1
    shell:
        op.as_one_line("""
        export CPLEX_DIR={params.CPLEX} && 
        JABBA_PATH=$(Rscript -e 'cat(paste0(installed.packages()["JaBbA", "LibPath"], "/JaBbA/extdata/"))'); $JABBA_PATH/jba `readlink -e {input.junc}` `readlink -e {input.rds}` --field foreground --cfield tier --outdir `dirname {output.rds}` --cores {threads} {params.opts} --gurobi TRUE > {log.stdout} 2> {log.stderr}
        """)


rule _jabba_link_graph_rds:
    input:
        rds = CFG["dirs"]["jabba"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/jabba.simple.gg.rds"
    output:
        rds = CFG["dirs"]["outputs"] + "rds/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.jabba.simple.gg.rds"
    wildcard_constraints:
        pair_status = "matched"
    run:
        op.relative_symlink(input.rds, output.rds)

# Generates the target sentinels for each run, which generate the symlinks
rule _jabba_all:
    input:
        expand(
            [
                str(rules._jabba_pon_link_dryclean_normal_rds.output.rds),
            ],
            zip,  # Run expand() with zip(), not product()
            seq_type=NORMALS["seq_type"],
            genome_build=NORMALS["genome_build"],
            sample_id=NORMALS["sample_id"],
            tissue_status=NORMALS["tissue_status"],
            allow_missing=True
            ),
        expand(
            [
                str(rules._jabba_link_graph_rds.output.rds),
            ],
            zip,
            seq_type=CFG["runs"]["tumour_seq_type"],
            genome_build=CFG["runs"]["tumour_genome_build"],
            tumour_id=CFG["runs"]["tumour_sample_id"],
            normal_id=CFG["runs"]["normal_sample_id"],
            pair_status=CFG["runs"]["pair_status"],
            allow_missing=True
        )


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)