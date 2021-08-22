#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Prasath Pararajalingam
# Contributors:     N/A


##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
import os.path

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
# `CFG` is a shortcut to `config["lcr-modules"]["sigprofiler"]`
CFG = op.setup_module(
    name = "sigprofiler",
    version = "1.0",
    subdirectories = ["inputs", "estimate", "extract", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _install_sigprofiler_genome,
    _sigprofiler_input_maf,
    _sigprofiler_output_tsv,
    _sigprofiler_all,


##### FUNCTIONS #####

def get_dir(wildcards):
    if "ID" in wildcards.type:
        topdir = 'ID'
    elif "SBS" in wildcards.type:
        topdir = 'SBS'
    elif "DBS" in wildcards.type:
        topdir = 'DBS'

    cc = config["lcr-modules"]["sigprofiler"]
    mat = cc["dirs"]["inputs"]+f"matrices/{wildcards.seq_type}--{wildcards.genome_build}/{wildcards.sample_set}/output/{topdir}/{wildcards.sample_set}.{wildcards.type}.all"
    ret = {'script' : cc["inputs"]["extractor"], 'mat' : mat}
    return(ret)

max_sigs = {
    "SBS6"    : 10,
    "SBS18"   : 10,
    "SBS24"   : 20,
    "SBS96"   : 20,
    "SBS288"  : 25,
    "SBS384"  : 25,
    "SBS4608" : 30,
    "SBS1536" : 30,
    "SBS6144" : 35,
    "ID28"    : 5,
    "ID83"    : 5,
    "ID96"    : 10,
    "ID332"   : 10,
    "ID415"   : 15,
    "ID8628"  : 20,
    "DBS78"   : 15,
    "DBS150"  : 15,
    "DBS186"  : 15,
    "DBS1248" : 20,
    "DBS2400" : 20,
    "DBS2976" : 20
}

##### RULES #####

rule _install_sigprofiler_genome:
    output:
        complete = CFG["dirs"]["inputs"] + "sigprofiler_prereqs/{genome_build}.installed"
    params:
        ref = lambda w: {"grch37":"GRCh37", "hg19":"GRCh37",
                         "grch38": "GRCh38", "hg38": "GRCh38"}[w.genome_build]
    conda: CFG["conda_envs"]["sigprofiler"]
    shell:
        op.as_one_line("""
        python -c 'from SigProfilerMatrixGenerator import install as genInstall;
        genInstall.install("{params.ref}", rsync = False, bash = True)'
            &&
        touch {output.complete}
        """)


# Symlinks the input files into the module results directory (under '00-inputs/')
rule _sigprofiler_input_maf:
    input:
        maf = CFG["inputs"]["maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "matrices/{seq_type}--{genome_build}/{sample_set}/{sample_set}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)

# Generates sample by k-mer context matrices from MAF
rule _sigprofiler_run_generator:
    input:
        mg = str(rules._install_sigprofiler_genome.output.complete),
        script = CFG["inputs"]["generator"],
        maf = str(rules._sigprofiler_input_maf.output.maf)
    output:
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/SBS/{{sample_set}}.{type}.all", type = ['SBS6','SBS18','SBS24','SBS96','SBS288','SBS384','SBS4608','SBS1536','SBS6144']),
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/ID/{{sample_set}}.{type}.all", type = ['ID28','ID83','ID96','ID332','ID415','ID8628']),
        expand(CFG["dirs"]["inputs"] + "matrices/{{seq_type}}--{{genome_build}}/{{sample_set}}/output/DBS/{{sample_set}}.{type}.all", type = ['DBS78','DBS150','DBS186','DBS1248','DBS2400','DBS2976'])
    params:
        ref = lambda w: {"grch37":"GRCh37", "hg19":"GRCh37",
                         "grch38": "GRCh38", "hg38": "GRCh38"}[w.genome_build]
    conda: CFG["conda_envs"]["sigprofiler"]
    threads: CFG["threads"]["generator"]
    shell:
        "python {input.script} {wildcards.sample_set} {params.ref} {input.maf}"

# Performs NMF on a wide range of ranks to obtain the optimal signature rank
rule _sigprofiler_run_estimate:
    input:
        unpack(get_dir)
    output:
        stat = CFG["dirs"]["estimate"]+"{seq_type}--{genome_build}/{sample_set}/{type}/All_solutions_stat.csv"
    params:
        ref = lambda w: {"grch37":"GRCh37", "hg19":"GRCh37",
                         "grch38": "GRCh38", "hg38": "GRCh38"}[w.genome_build],
        context_type = '96,DINUC,ID',
        exome = lambda w: {'genome': 'False', 'capture': 'True'}[w.seq_type],
        min_sig = 1,
        max_sig = lambda w: max_sigs[w.type],
        nmf_repl = 30,
        norm = 'gmm',
        nmf_init = 'nndsvd_min',
        outpath = CFG["dirs"]["estimate"]+"{seq_type}--{genome_build}/{sample_set}"
    conda: CFG["conda_envs"]["sigprofiler"]
    threads: CFG["threads"]["estimate"]
    shell:
        op.as_one_line("""
        python {input.script} {input.mat} {params.outpath} 
        {params.ref} {params.context_type} {params.exome} {params.min_sig} {params.max_sig} 
        {params.nmf_repl} {params.norm} {params.nmf_init} {threads}
        """)

# Runs NMF around the optimal rank from _sigprofiler_run_estimate
# using higher NMF replications to get final signatures
rule _sigprofiler_run_extract:
    input:
        unpack(get_dir),
        stat = str(rules._sigprofiler_run_estimate.output.stat)
    output:
        decomp = CFG["dirs"]["extract"]+"{seq_type}--{genome_build}/{sample_set}/{type}/Suggested_Solution/COSMIC_{type}_Decomposed_Solution/De_Novo_map_to_COSMIC_{type}.csv"
    params:
        ref = lambda w: {"grch37":"GRCh37", "hg19":"GRCh37", 
                         "grch38": "GRCh38", "hg38": "GRCh38"}[w.genome_build],
        context_type = '96,DINUC,ID',
        exome = lambda w: {'genome': 'False', 'capture': 'True'}[w.seq_type],
        nmf_repl = 200,
        norm = 'gmm',
        nmf_init = 'nndsvd_min',
        outpath = CFG["dirs"]["extract"]+"{seq_type}--{genome_build}/{sample_set}"
    conda: CFG["conda_envs"]["sigprofiler"]
    threads: CFG["threads"]["extract"]
    shell:
        op.as_one_line("""
        python {input.script} {input.mat} {params.outpath} 
        {params.ref} {params.context_type} {params.exome}
        $(awk 'BEGIN {{OFS=FS=","}} $1 ~ "*" {{S=substr($1,1,length($1)-1); if (S-2<1) {{print 1}} else {{print S-2}}}}' {input.stat})
        $(awk 'BEGIN {{OFS=FS=","}} $1 ~ "*" {{S=substr($1,1,length($1)-1); print S+2}}' {input.stat})
        {params.nmf_repl} {params.norm} {params.nmf_init} {threads}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _sigprofiler_output_tsv:
    input:
        decomp = str(rules._sigprofiler_run_extract.output.decomp)
    output:
        decomp = CFG["dirs"]["outputs"] + "cosmic_sigs/{seq_type}--{genome_build}/{sample_set}.COSMIC_{type}.De_Novo_map_to_COSMIC_{type}.csv"
    run:
        op.relative_symlink(input.decomp, output.decomp)

# Generates the target sentinels for each run, which generate the symlinks
rule _sigprofiler_all:
    input:
        expand(
            expand(
                CFG["dirs"]["outputs"] + "cosmic_sigs/{seq_type}--{genome_build}/{sample_set}.COSMIC_{{type}}.De_Novo_map_to_COSMIC_{{type}}.csv",
                zip,
                seq_type = CFG["mafs"]["seq_type"],
                genome_build = CFG["mafs"]["projected_build"],
                sample_set = CFG["mafs"]["sample_set"]),
            type = CFG["inputs"]["type"])


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
