#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Cancer IT
# Module Author:    Helena Winata
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
# `CFG` is a shortcut to `config["lcr-modules"]["cluster_sv"]`
CFG = op.setup_module(
    name = "cluster_sv",
    version = "1.0",
    subdirectories = ["inputs", "reformat_bedpe", "cluster_sv", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _cluster_sv_input_bedpe,
    _cluster_sv_reformat_bedpe,
    _cluster_sv_output_tsv,
    _cluster_sv_all,


##### RULES #####

# download cluster files from git repo without cloning the repo itself
# decompress files into the 00-inputs
rule _cluster_sv_install:
    output:
        cluster_sv = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/R/run_cluster_sv.R", # one of the files in the repo
        chr_sizes = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/references/{genome_build}.chrom_sizes", # another script from the repo
        coords = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/references/{genome_build}_centromere_and_telomere_coords.txt" # convert output to bed format
    params:
        url = "https://github.com/whelena/ClusterSV/archive/refs/tags/v" + str(CFG["options"]["cluster_sv_version"]) + ".tar.gz",
        folder = CFG["dirs"]["inputs"]
    shell:
        op.as_one_line("""
        wget -qO- {params.url} |
        tar xzf - -C {params.folder}
        """)

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _cluster_sv_input_bedpe:
    input:
        bedpe = CFG["inputs"]["sample_bedpe"]
    output:
        bedpe = CFG["dirs"]["inputs"] + "bedpe/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    run:
        op.relative_symlink(input.bedpe, output.bedpe)


# Example variant calling rule (multi-threaded; must be run on compute server/cluster)
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

''''
def _cluster_sv_get_references(wildcards):
    CFG = config["lcr-modules"]["cluster_sv"]
    if CFG["options"]["chr"] is not None:
        ref1 = CFG["options"]["chr"]
    elif CFG["options"]["chr"] is None:
        ref1 = chr_sizes = CFG["dirs"]["inputs"] + "ClusterSV-" + str(CFG["options"]["cluster_sv_version"]) + "/references/{genome_build}.chrom_sizes"



        def _cluster_sv_get_references_custom(wildcards, ref):
            if CFG["options"][ref] is None:
'''


rule _cluster_sv_run:
    input:
        bedpe = str(rules._cluster_sv_reformat_bedpe.output.bedpe),
        cluster_sv = str(rules._cluster_sv_install.output.cluster_sv),
        chr_sizes = str(rules._cluster_sv_install.output.chr_sizes),
        coords = str(rules._cluster_sv_install.output.coords)
    output:
        outdir = directory(CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}"),
        clusters = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_clusters_and_footprints.tsv",
        pval = CFG["dirs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/sv_distance_pvals.tsv"
    log:
        stdout = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stdout.log",
        stderr = CFG["logs"]["cluster_sv"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/cluster_sv.stderr.log"
    # params:
    #     cluster_sv_script = CFG["scripts"]["cluster_sv"],
    #     chr_sizes = CFG["options"]["chr_sizes"],
    #     coords = CFG["options"]["cen_telo_coords"]
    conda:
        CFG["conda_envs"]["cluster_sv"]
    threads:
        CFG["threads"]["cluster_sv"]
    resources:
        mem_mb = CFG["mem_mb"]["cluster_sv"]
    shell:
        op.as_one_line("""
        Rscript {input.cluster_sv} 
        -bedpe {input.bedpe} 
        -chr {input.chr_sizes} 
        -cen_telo {input.coords} 
        -out {output.outdir} 
        -n {threads}
        > {log.stdout} 2> {log.stderr} &&
        sed -i '1i #CHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tID\tQUAL\tSTRAND_A\tSTRAND_B\tCLUSTER_ID\tNUM_SV\tFP_ID_LOW\tFP_ID_HIGH\tFP_COORDS_LOW\tFP_COORDS_HIGH\tP_VAL' {output.clusters}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _cluster_sv_output_tsv:
    input:
        clusters = str(rules._cluster_sv_run.output.pval),
        pval = str(rules._cluster_sv_run.output.pval)
    output:
        clusters = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_clusters_and_footprints.tsv",
        pval = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.sv_distance_pvals.tsv"
    run:
        op.relative_symlink(input.clusters, output.clusters)
        op.relative_symlink(input.pval, output.pval)


# Generates the target sentinels for each run, which generate the symlinks
rule _cluster_sv_all:
    input:
        expand(
            [
                str(rules._cluster_sv_output_tsv.output.clusters),
                str(rules._cluster_sv_output_tsv.output.pval)
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
