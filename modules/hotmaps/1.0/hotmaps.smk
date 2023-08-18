#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
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
# `CFG` is a shortcut to `config["lcr-modules"]["hotmaps"]`
CFG = op.setup_module(
    name = "hotmaps",
    version = "1.0",
    subdirectories = ["inputs", "hotmaps", "outputs"],
)

# Define rules to be run locally when using a compute cluster
# TODO: Replace with actual rules once you change the rule names
localrules:
    _hotmaps_input_maf,
    _hotmaps_step_2,
    _hotmaps_output_txt,
    _hotmaps_all,


##### RULES #####

rule _install_hotmaps:
    output:
        complete = CFG["dirs"]["inputs"] + "hotmaps_installed.complete"
    params:
        hotmaps_repo = CFG["hotmaps_repo"],
        output_dir = CFG["dirs"]["inputs"]
    shell:
        op.as_one_line("""
        wget -qcP {params.output_dir} {params.hotmaps_repo} &&
        unzip -d {params.output_dir} {params.output_dir}/master.zip &&
        rm {params.output_dir}/master.zip &&
        touch {output.complete}
        """)

rule _hotmaps_update_config:
    input:
        hotmaps_installed = str(rules._install_hotmaps.output.complete)
    output:
        hotmaps_config_updated = CFG["dirs"]["inputs"] + "config_updated.complete"
    params:
        modbase_dir = "modbase_dir=" + CFG["modbase_dir"],
        pdb_dir = "pdb_dir=" + CFG["pdb_dir"],
        refseq_homology_dir = "refseq_homology: " + CFG["refseq_homology_dir"],
        ensembl_homology_dir = "ensembl_homology: " + CFG["ensembl_homology_dir"],
        biological_assembly_dir = "biological_assembly: " + CFG["biological_assembly_dir"],
        non_biological_assembly_dir = "non_biological_assembly: " + CFG["non_biological_assembly_dir"],
        config_file = CFG["dirs"]["inputs"] + "HotMAPS-master/config.txt"
    shell:
        op.as_one_line("""
        sed -E -i "s@(modbase_dir=).*@{params.modbase_dir}@" {params.config_file} &&
        sed -E -i "s@(pdb_dir=).*@{params.pdb_dir}@" {params.config_file} &&
        sed -E -i "s@(refseq_homology:).*@{params.refseq_homology_dir}@" {params.config_file} &&
        sed -E -i "s@(ensembl_homology:).*@{params.ensembl_homology_dir}@" {params.config_file} &&
        sed -E -i "s@(biological_assembly:).*@{params.biological_assembly_dir}@" {params.config_file} &&
        sed -E -i "s@(non_biological_assembly:).*@{params.non_biological_assembly_dir}@" {params.config_file} &&
        touch {output.hotmaps_config_updated}
        """)

rule _hotmaps_get_pdb_info:
    input:
        hotmaps_installed = str(rules._install_hotmaps.output.complete)
    output:
        pdb_info = CFG["dirs"]["inputs"] + "pdb/pdb_info.txt"
    params:
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_passwd = CFG["options"]["mysql"]["mysql_passwd"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_database = CFG["options"]["mysql"]["mysql_db"],
        mysql_db_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/sql/get_pdb_info.sql"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        mysql -u {params.mysql_user} -A -p{params.mysql_passwd} -h {params.mysql_host} {params.mysql_database} <
        {params.mysql_db_script} > {output.pdb_info}
        """)

rule _hotmaps_add_pdb_path:
    input:
        pdb_info = str(rules._hotmaps_get_pdb_info.output.pdb_info),
        config = str(rules._hotmaps_update_config.output.hotmaps_config_updated)
    output:
        pdb_path = CFG["dirs"]["inputs"] + "pdb/pdb_info.path.txt"
    params:
        add_path_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/add_path_info.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        "python {params.add_path_script} -p {input.pdb_info} -o {output.pdb_path}"

rule _hotmaps_add_pdb_description:
    input:
        pdb_info_path = str(rules._hotmaps_add_pdb_path.output.pdb_path)
    output:
        fully_described_pdb = CFG["dirs"]["inputs"] + "pdb/fully_described_pdb_info.txt"
    params:
        describe_pdb_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/chain_description.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        "python {params.describe_pdb_script} -i {input.pdb_info_path} -o {output.fully_described_pdb}"

# Symlinks the input files into the module results directory (under '00-inputs/')
#rule _hotmaps_input_maf:
#    input:
#        maf = CFG["inputs"]["sample_maf"]
#    output:
#        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{sample_id}.maf"
#    group: 
#        "input_and_step_1"
#    run:
#        op.absolute_symlink(input.maf, output.maf)


## Example variant calling rule (multi-threaded; must be run on compute server/#cluster)
## TODO: Replace example rule below with actual rule
#rule _hotmaps_step_1:
#    input:
#        tumour_maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}.maf",
#        normal_maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{normal_id}.maf",
#        fasta = reference_files("genomes/{genome_build}/genome_fasta/genome.fa")
#    output:
#        txt = CFG["dirs"]["hotmaps"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.txt"
#    log:
#        stdout = CFG["logs"]["hotmaps"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stdout.log",
#        stderr = CFG["logs"]["hotmaps"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_1.stderr.log"
#    params:
#        opts = CFG["options"]["step_1"]
#    conda:
#        CFG["conda_envs"]["samtools"]
#    threads:
#        CFG["threads"]["step_1"]
#    resources:
#        **CFG["resources"]["step_1"]
#    group: 
#        "input_and_step_1"
#    shell:
#        op.as_one_line("""
#        <TODO> {params.opts} --tumour {input.tumour_maf} --normal {input.normal_maf}
#        --ref-fasta {input.fasta} --output {output.txt} --threads {threads}
#        > {log.stdout} 2> {log.stderr}
#        """)


# Example variant filtering rule (single-threaded; can be run on cluster head node)
# TODO: Replace example rule below with actual rule
#rule _hotmaps_step_2:
#    input:
#        txt = str(rules._hotmaps_step_1.output.txt)
#    output:
#        txt = CFG["dirs"]["hotmaps"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/output.filt.txt"
#    log:
#        stderr = CFG["logs"]["hotmaps"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/step_2.stderr.log"
#    params:
#        opts = CFG["options"]["step_2"]
#    shell:
#        "grep {params.opts} {input.txt} > {output.txt} 2> {log.stderr}"


## Symlinks the final output files into the module results directory (under #'99-outputs/')
## TODO: If applicable, add an output rule for each file meant to be exposed to #the user
#rule _hotmaps_output_txt:
#    input:
#        txt = str(rules._hotmaps_step_2.output.txt)
#    output:
#        txt = CFG["dirs"]["outputs"] + "txt/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.output.filt.txt"
#    run:
#        op.relative_symlink(input.txt, output.txt, in_module= True)


# Generates the target sentinels for each run, which generate the symlinks
rule _hotmaps_all:
    input:
        str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
