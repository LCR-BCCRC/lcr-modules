#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Cancer IT
# Module Author:    Helena Winata
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["cluster_sv"]`
CFG = op.setup_module(
    name = "cluster_sv",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "reformat_bedpe", "cluster_sv", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _cluster_sv_input_bedpe,
    _cluster_sv_reformat_bedpe,
    _cluster_sv_output_bedpe,
    _cluster_sv_all,


##### RULES #####


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cluster_sv_input_bedpe:
    input:
        bedpe = CFG["inputs"]["sample_bedpe"]
    output:
        bedpe = CFG["dirs"]["inputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
# TODO: Replace example rule below with actual rule
rule _cluster_sv_reformat_bedpe:
    input:
        bedpe = str(rules._cluster_sv_input_bedpe.output.bedpe)
    output:
        bedpe = CFG["dirs"]["reformat_bedpe"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    shell:
        op.as_one_line("""
        sed '/^#/d' {input.bedpe} | 
        sed -e 's/chr//g' |
        awk -F '\t' 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' 
        > {output.bedpe}
        """)


rule _cluster_sv_run:
    input:
        bedpe = str(rules._cluster_sv_reformat_bedpe.output.bedpe)
        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
    output:
        outdir = directory(CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}"),
        clusters = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_clusters_and_footprints.tsv",
        pval = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_distance_pvals.tsv"
    log:
        stdout = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stdout.log",
        stderr = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stderr.log"
    params:
        cluster_sv_script = CFG["scripts"]["cluster_sv"]
    conda:
        CFG["conda_envs"]["cluster_sv"]
    threads:
        CFG["threads"]["cluster_sv"]
    resources:
        mem_mb = CFG["mem_mb"]["cluster_sv"]
    shell: #TODO
        op.as_one_line("""
        Rscript {params.cluster_sv_script} -n {threads}
        > {log.stdout} 2> {log.stderr}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
# TODO: If applicable, add an output rule for each file meant to be exposed to the user
rule _cluster_sv_output_bedpe:
    input:
        bedpe = str(rules._cluster_sv_run.output.bedpe)
    output:
        bedpe = CFG["dirs"]["outputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)


# Generates the target sentinels for each run, which generate the symlinks
rule _cluster_sv_all:
    input:
        expand(
            [
                str(rules._cluster_sv_output_bedpe.output.bedpe),
                # TODO: If applicable, add other output rules here
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
