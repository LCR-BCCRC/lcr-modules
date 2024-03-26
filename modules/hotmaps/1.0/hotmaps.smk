#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


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
# `CFG` is a shortcut to `config["lcr-modules"]["hotmaps"]`
CFG = op.setup_module(
    name = "hotmaps",
    version = "1.0",
    subdirectories = ["inputs", "maf2vcf", "bcftools", "vcf2maf", "hotmaps", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _hotmaps_input_maf,
    _hotmaps_input_subsetting_categories,
    _hotmaps_prep_input,
    _hotmaps_split_dnps,
    _hotmaps_maf2vcf,
    _hotmaps_bcftools,
    _hotmaps_vcf2maf,
    _hotmaps_merge_mafs,
    _hotmaps_deblacklist,
    _hotmaps_input,
    _install_hotmaps,
    _hotmaps_update_config,
    _hotmaps_get_pdb_info,
    _hotmaps_add_pdb_path,
    _hotmaps_add_pdb_description,
    _hotmaps_prep_mutations,
    _hotmaps_prep_mupit_annotation,
    _hotmaps_filter_hypermutated,
    _hotmaps_count_mutations,
    _hotmaps_format_mutations,
    _hotmaps_merge_mutations,
    _hotmaps_get_mutations,
    _hotmaps_split_pdbs,
    _hotmaps_merge_hotspots,
    _hotmaps_multiple_test_correct,
    _hotmaps_find_gene,
    _hotmaps_find_structure,
    _hotmaps_output,
    _hotmaps_detailed_hotspots,
    _hotmaps_aggregate,
    _hotmaps_all,

##### RULES #####

if "launch_date" in CFG:
    launch_date = CFG["launch_date"]
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS = os.path.abspath(config["lcr-modules"]["hotmaps"]["maf_processing"]["prepare_mafs"])

# Set path to vcf2maf script using config value
VCF2MAF_SCRIPT_PATH = CFG["options"]["vcf2maf"]["src_dir"]

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _hotmaps_input_maf:
    input:
        maf = CFG["inputs"]["input_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "master_maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

rule _hotmaps_input_subsetting_categories:
    input:
        subsetting_categories = ancient(CFG["inputs"]["subsetting_categories"])
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

checkpoint _hotmaps_prep_input:
    input:
        maf = expand(
            str(rules._hotmaps_input_maf.output.maf),
            seq_type = CFG["samples"]["seq_type"].unique(),
            allow_missing = True
            ),
        subsetting_categories = ancient(str(rules._hotmaps_input_subsetting_categories.output.subsetting_categories))
    output:
        CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/done"
    log:
        stdout = CFG["logs"]["inputs"] + "{genome_build}/{sample_set}--{launch_date}/prep_input/prep_input_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["maf_processing"]["include_non_coding"]).upper(),
        mode = "HotMAPS",
        metadata_cols = CFG["samples"],
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

rule _hotmaps_split_dnps:
    input:
        maf = CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.maf"
    output:
        dnps = temp(CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.dnps.maf"),
        completed = CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}_split_dnps.complete"
    shell:
        op.as_one_line("""
        variant_type_col=$(head -n 1 {input.maf} | sed 's/\\t/\\n/g' | nl | grep "Variant_Type" | cut -f 1) &&
        protein_position_col=$(head -n 1 {input.maf} | sed 's/\\t/\\n/g' | nl | grep "Protein_position" | cut -f 1) &&
        cat <( head -n 1 {input.maf} ) <( awk -v var_col="$variant_type_col" -v protein_col="$protein_position_col" ' {{ if ( $var_col=="DNP" && $protein_col ~ /[0-9?]+[-][0-9?]+/) print $0 }} ' {input.maf} ) > {output.dnps} &&
        touch {output.completed}
        """)

checkpoint _hotmaps_maf2vcf:
    input:
        dnps = str(rules._hotmaps_split_dnps.output.completed),
        fasta = reference_files("genomes/grch37/genome_fasta/genome.fa")
    output:
        done = CFG["dirs"]["maf2vcf"] + "completed/{genome_build}/{sample_set}--{launch_date}/{md5sum}.done"
    params:
        vcf = CFG["dirs"]["maf2vcf"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.dnps.vcf",
        vcf_dir = CFG["dirs"]["maf2vcf"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/",
        dnps = str(rules._hotmaps_split_dnps.output.dnps)
    conda:
        CFG["conda_envs"]["vcf2maf"]
    log:
        stdout = CFG["logs"]["maf2vcf"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}_maf2vcf.stdout.log",
        stderr = CFG["logs"]["maf2vcf"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}_maf2vcf.stderr.log"
    shell:
        op.as_one_line("""
        if [ $( wc -l < {params.dnps}) -gt 1 ] ;
        then
            mkdir -p {params.vcf_dir} &&
            maf2vcf.pl 
            --input-maf {params.dnps} 
            --output-dir {params.vcf_dir} 
            --output-vcf {params.vcf} 
            --ref-fasta {input.fasta} 
            --per-tn-vcfs 
            > {log.stdout} 2> {log.stderr} ;
        fi &&
        touch {output.done}
        """)

rule _hotmaps_bcftools:
    input:
        vcf = CFG["dirs"]["maf2vcf"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}.vcf",
        dnp_vcf = str(rules._hotmaps_maf2vcf.output.done)
    output:
        vcf = CFG["dirs"]["bcftools"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}.annotate.vcf"
    conda:
        CFG["conda_envs"]["bcftools"] 
    log:
        stdout = CFG["logs"]["bcftools"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}/bcftools_norm.stdout.log",
        stderr = CFG["logs"]["bcftools"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}/bcftools_norm.stderr.log"
    shell:
        op.as_one_line("""
        bcftools norm --atomize 
        --output {output.vcf} 
        --output-type v {input.vcf} 
        > {log.stdout} 2> {log.stderr}
        """)

VCF2MAF_VERSION_MAP = {
    "grch37":"GRCh37"
}

rule _hotmaps_vcf2maf:
    input:
        vcf = str(rules._hotmaps_bcftools.output.vcf),
        fasta = reference_files("genomes/grch37/genome_fasta/genome.fa"),
        vep_cache = CFG["options"]["vcf2maf"]["vep_cache"]
    output:
        maf = CFG["dirs"]["vcf2maf"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}.maf",
        vep = CFG["dirs"]["bcftools"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}.annotate.vep.vcf"
    params:
        opts = CFG["options"]["vcf2maf"]["options"],
        build = lambda w: VCF2MAF_VERSION_MAP[w.genome_build],
        custom_enst = lambda w: "--custom-enst " + str(config["lcr-modules"]["hotmaps"]["options"]["vcf2maf"]["custom_enst"][w.genome_build]) if config["lcr-modules"]["hotmaps"]["options"]["vcf2maf"]["custom_enst"][w.genome_build] is not None else ""
    conda:
        CFG["conda_envs"]["vcf2maf"]
    log:
        stdout = CFG["logs"]["vcf2maf"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}/vcf2maf.stdout.log",
        stderr = CFG["logs"]["vcf2maf"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/{tumour_id}_vs_{normal_sample_id}/vcf2maf.stderr.log"
    shell:
        op.as_one_line("""
        VCF2MAF_SCRIPT_PATH={VCF2MAF_SCRIPT_PATH};
        PATH=$VCF2MAF_SCRIPT_PATH:$PATH;
        VCF2MAF_SCRIPT="$VCF2MAF_SCRIPT_PATH/vcf2maf.pl";
        vepPATH=$(dirname $(which variant_effect_predictor.pl))/../share/variant-effect-predictor* ;
        if [[ $(which vcf2maf.pl) =~ $(ls $VCF2MAF_SCRIPT) ]]; then
            echo "using bundled patched script $VCF2MAF_SCRIPT";
            echo "Using $VCF2MAF_SCRIPT to run {rule} for {wildcards.tumour_id} on $(hostname) at $(date)" > {log.stderr};
            vcf2maf.pl 
            --input-vcf {input.vcf} 
            --output-maf {output.maf} 
            --tumor-id {wildcards.tumour_id} 
            --normal-id {wildcards.normal_sample_id} 
            --vcf-tumor-id TUMOR 
            --vcf-normal-id NORMAL 
            --ncbi-build {params.build} 
            --vep-path $vepPATH 
            --vep-data {input.vep_cache}
            --ref-fasta {input.fasta}
            --retain-info gnomADg_AF,gnomADg_AF
            {params.opts} {params.custom_enst}
            >> {log.stdout} 2>> {log.stderr};
            if [[ $(wc -l < {output.maf}) -eq 0 ]]; then
                echo "vcf2maf failed for {wildcards.tumour_id}, running again" >> {log.stderr};
                vcf2maf.pl
                --input-vcf {input.vcf} 
                --output-maf {output.maf}
                --tumor-id {wildcards.tumour_id} 
                --normal-id {wildcards.normal_sample_id} 
                --vcf-tumor-id TUMOR 
                --vcf-normal-id NORMAL 
                --ncbi-build {params.build} 
                --vep-path $vepPATH 
                --ref-fasta {input.fasta} 
                --retain-info gnomADg_AF,gnomADg_AF
                {params.opts} {params.custom_enst}
                >> {log.stdout} 2>> {log.stderr};
            fi;
            if [[ $(wc -l < {output.maf}) -eq 0 ]]; then
                echo "Could not resolve error running {wildcards.tumour_id}" >> {log.stderr};
                exit 1;
            fi;
        else echo "ERROR: PATH is not set properly, using $(which vcf2maf.pl) will result in error during execution. Please ensure $VCF2MAF_SCRIPT exists." > {log.stderr}; fi
        """) 

def _get_dnp_mafs(wildcards):
    CFG = config["lcr-modules"]["hotmaps"]
    checkpoint_outputs = checkpoints._hotmaps_maf2vcf.get(**wildcards).output.done

    maf2vcf_output_dir = expand(
        CFG["dirs"]["maf2vcf"] + "vcf/{genome_build}/{sample_set}--{launch_date}/{md5sum}/", 
        genome_build = "grch37",
        sample_set = wildcards.sample_set,
        launch_date = wildcards.launch_date,
        md5sum = wildcards.md5sum)
    
    TUMOUR_IDS, NORMAL_IDS = glob_wildcards(maf2vcf_output_dir[0] + "{tumour_id}_vs_{normal_sample_id}.vcf")

    return expand(
        str(rules._hotmaps_vcf2maf.output.maf),
        zip,
        tumour_id = TUMOUR_IDS,
        normal_sample_id = NORMAL_IDS,
        allow_missing = True
    )

rule _hotmaps_merge_mafs:
    input:
        maf_annotated = _get_dnp_mafs,
        maf2vcf = str(rules._hotmaps_maf2vcf.output.done),
        maf_created = str(rules._hotmaps_prep_input.output)
    output:
        maf = temp(CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.reannotated.maf")
    params:
        maf = CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.maf"
    log:
        stderr = CFG["logs"]["inputs"] + "merge/{genome_build}/{sample_set}--{launch_date}/{md5sum}.stderr.log"
    run:
        import pandas as pd
        error_out = open(log.stderr, "w")
        errors = 0
        main_maf = pd.read_table(params.maf, comment = "#", sep="\t")
        main_maf = main_maf[~((main_maf["Variant_Type"]=="DNP") & (main_maf["Protein_position"].str.match(r'\d+-\d+|\?-\d+|\d+-\?')))]
        for maf in input.maf_annotated:
            try:
                df = pd.read_table(maf, comment="#", sep="\t")
                main_maf = pd.concat([main_maf, df])
            except:
                error_out.write(f"Error reading MAF: {maf}\n")
                errors += 1
        if errors == 0:
            main_maf.to_csv(output.maf, sep="\t", na_rep="NA", index=False)
        error_out.close()

rule _hotmaps_deblacklist:
    input:
        maf = ancient(str(rules._hotmaps_merge_mafs.output.maf)),
        dnps = str(rules._hotmaps_split_dnps.output.dnps),
        blacklists = CFG["maf_processing"]["blacklists"],
        deblacklist_script = CFG["deblacklist_script"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/{md5sum}.reannotated.deblacklisted.maf"
    params:
        drop_threshold = CFG["maf_processing"]["blacklist_drop_threshold"],
        blacklists = CFG["maf_processing"]["blacklists"]
    log:
        stdout = CFG["logs"]["inputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stdout.log",
        stderr = CFG["logs"]["inputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/deblacklist/deblacklist.stderr.log"
    shell:
        op.as_one_line("""
        {input.deblacklist_script} 
        --input {input.maf} 
        --output {output.maf} 
        --drop-threshold {params.drop_threshold} 
        --blacklists {params.blacklists} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_input:
    input:
        maf = str(rules._hotmaps_deblacklist.output.maf)
    output:
        maf = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/input.{sample_set}.maf"
    run:
        import pandas as pd
        maf = pd.read_table(input.maf, comment="#", sep="\t")
        maf = maf[["Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Hugo_Symbol","Variant_Classification","HGVSp_Short","Transcript_ID","Strand"]]
        maf.to_csv(output.maf, sep="\t", na_rep="NA", index=False)

rule _install_hotmaps:
    output:
        installed = CFG["dirs"]["inputs"] + "HotMAPS-master/hotmaps_installed.done",
        config = CFG["dirs"]["inputs"] + "HotMAPS-master/config.txt"
    params:
        hotmaps_repo = CFG["hotmaps_repo"],
        output_dir = CFG["dirs"]["inputs"]
    shell:
        op.as_one_line("""
        wget -qcP {params.output_dir} {params.hotmaps_repo} &&
        unzip -d {params.output_dir} {params.output_dir}/master.zip &&
        rm {params.output_dir}/master.zip &&
        touch {output.installed}
        """)

# This updates the config file read by hotmaps with the paths pointing to PDB structure locations
rule _hotmaps_update_config:
    input:
        hotmaps_installed = str(rules._install_hotmaps.output.installed),
        config = str(rules._install_hotmaps.output.config)
    output:
        config_updated = CFG["dirs"]["inputs"] + "HotMAPS-master/config_update.done"
    params:
        modbase_dir = "modbase_dir=" + CFG["pdb_structure_dirs"]["modbase_dir"],
        pdb_dir = "pdb_dir=" + CFG["pdb_structure_dirs"]["pdb_dir"],
        refseq_homology_dir = "refseq_homology: " + CFG["pdb_structure_dirs"]["refseq_homology_dir"],
        ensembl_homology_dir = "ensembl_homology: " + CFG["pdb_structure_dirs"]["ensembl_homology_dir"],
        biological_assembly_dir = "biological_assembly: " + CFG["pdb_structure_dirs"]["biological_assembly_dir"],
        non_biological_assembly_dir = "non_biological_assembly: " + CFG["pdb_structure_dirs"]["non_biological_assembly_dir"]
    shell:
        op.as_one_line("""
        sed -E -i "s@(modbase_dir=).*@{params.modbase_dir}@" {input.config} &&
        sed -E -i "s@(pdb_dir=).*@{params.pdb_dir}@" {input.config} &&
        sed -E -i "s@(refseq_homology:).*@{params.refseq_homology_dir}@" {input.config} &&
        sed -E -i "s@(ensembl_homology:).*@{params.ensembl_homology_dir}@" {input.config} &&
        sed -E -i "s@(biological_assembly:).*@{params.biological_assembly_dir}@" {input.config} &&
        sed -E -i "s@(non_biological_assembly:).*@{params.non_biological_assembly_dir}@" {input.config} &&
        touch {output.config_updated}
        """)

rule _hotmaps_get_pdb_info:
    input:
        hotmaps_installed = str(rules._install_hotmaps.output.installed)
    output:
        pdb_info = temp(CFG["dirs"]["inputs"] + "pdb/pdb_info.txt")
    params:
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_passwd = CFG["options"]["mysql"]["mysql_passwd"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_database = CFG["options"]["mysql"]["mysql_db"],
        mysql_db_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/sql/get_pdb_info.sql"
    log:
        stderr = CFG["logs"]["inputs"] + "get_pdb_info/get_pdb_info.stderr.log"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        mysql -u {params.mysql_user} 
        -A -p{params.mysql_passwd} 
        -h {params.mysql_host} {params.mysql_database} 
        < {params.mysql_db_script} > {output.pdb_info} 
        2> {log.stderr}
        """)

rule _hotmaps_add_pdb_path:
    input:
        pdb_info = str(rules._hotmaps_get_pdb_info.output.pdb_info),
        config = str(rules._hotmaps_update_config.output.config_updated)
    output:
        pdb_path = temp(CFG["dirs"]["inputs"] + "pdb/pdb_info.path.txt")
    params:
        add_path_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/add_path_info.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["inputs"] + "add_pdb_path/add_pdb_path.stdout.log",
        stderr = CFG["logs"]["inputs"] + "add_pdb_path/add_pdb_path.stderr.log"
    shell:
        op.as_one_line("""
        python {params.add_path_script} 
        -p {input.pdb_info} 
        -o {output.pdb_path} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_add_pdb_description:
    input:
        pdb_info_path = str(rules._hotmaps_add_pdb_path.output.pdb_path)
    output:
        fully_described_pdb = CFG["dirs"]["inputs"] + "pdb/fully_described_pdb_info.txt"
    params:
        describe_pdb_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/chain_description.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["inputs"] + "add_pdb_description/add_pdb_description.stdout.log",
        stderr = CFG["logs"]["inputs"] + "add_pdb_description/add_pdb_description.stderr.log"
    shell:
        op.as_one_line("""
        python {params.describe_pdb_script} 
        -i {input.pdb_info_path} 
        -o {output.fully_described_pdb} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_prep_mutations:
    input:
        maf = str(rules._hotmaps_input.output.maf),
        installed = str(rules._install_hotmaps.output.installed)
    output:
        mupit_non_filtered = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/non_filtered_mupit.input.{sample_set}.maf"
    params:
        mut_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/",
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"],
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mut_regex = "'^input.+\.maf'$",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/map_maf_to_structure.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/prep_mutations/prep_mutations.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/prep_mutations/prep_mutations.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        --data-dir {params.mut_dir} 
        --match-regex {params.mut_regex} 
        --host {params.mysql_host} 
        --db {params.mysql_db} 
        --mysql-user {params.mysql_user} 
        --mysql-passwd {params.mysql_pass} 
        --output-dir {params.mut_dir} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_prep_mupit_annotation:
    input:
        maf = str(rules._hotmaps_input.output.maf)
    output:
        annotation = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mupit_annotations/mupit_mutations_{sample_set}"
    params:
        mut_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/",
        tumor_type = "{sample_set}",
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"],
        cov_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/",
        hypermut = "-mt " + CFG["options"]["hotmaps"]["hypermut_threshold"] if CFG["options"]["hotmaps"]["hypermut_threshold"] is not None else "",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/maf/convert_maf_to_mupit.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/prep_mupit_annotation/prep_mupit_annotation.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/prep_mupit_annotation/prep_mupit_annotation.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        --maf {input.maf} 
        -mh {params.mysql_host} 
        -mdb {params.mysql_db} 
        --mysql-user {params.mysql_user} 
        --mysql-passwd {params.mysql_pass} 
        --tumor-type {params.tumor_type} 
        --no-stratify {params.hypermut} 
        -i {params.cov_dir} 
        --output {output.annotation} 
        > {log.stdout} 2> {log.stderr}
        """)

# Removes hypermutated samples from analysis to minimize effects of passenger mutations based on Kandoth et al. (2013) paper
rule _hotmaps_filter_hypermutated:
    input:
        maf = str(rules._hotmaps_input.output.maf),
        nf_mupit = str(rules._hotmaps_prep_mutations.output.mupit_non_filtered)
    output:
        hypermut = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/hypermutated.{sample_set}.txt",
        mupit = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/mupit.input.{sample_set}.maf"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/filter_hypermutated.py",
        mut_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/",
        mut_regex = "'^input.+\.maf'$",
        mut_threshold = "--mut-threshold " + CFG["options"]["hotmaps"]["hypermut_threshold"] if CFG["options"]["hotmaps"]["hypermut_threshold"] is not None else "",
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/filter_hypermutated/filter_hypermutated.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/filter_hypermutated/filter_hypermutated.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        --raw-dir {params.mut_dir} 
        --match-regex {params.mut_regex} 
        {params.mut_threshold} 
        --sample-col Tumor_Sample_Barcode 
        --data-dir {params.mut_dir} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_count_mutations:
    input:
        mupit = str(rules._hotmaps_filter_hypermutated.output.mupit)
    output:
        collected = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/collected.input.{sample_set}.maf"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/count_mutations.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/"
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/count_mutations/count_mutations.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/count_mutations/count_mutations.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} --data-dir {params.data_dir} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_format_mutations:
    input:
        collected = str(rules._hotmaps_count_mutations.output.collected)
    output:
        tcga = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/mutation_tcga.{sample_set}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/format_mutations_table.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/"
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/format_mutations/format_mutations.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/format_mutations/format_mutations.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} --data-dir {params.data_dir} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_merge_mutations:
    input:
        tcga = str(rules._hotmaps_format_mutations.output.tcga)
    output:
        mysql_mut = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/mysql.mutations.tcga.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/merge_mutations_table_data.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/"
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/merge_mutations/merge_mutations.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/merge_mutations/merge_mutations.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} {params.data_dir} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_get_mutations:
    input:
        mysql_mut = str(rules._hotmaps_merge_mutations.output.mysql_mut)
    output:
        mutations = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mutations/mutations.txt"
    params:
        load_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/load_mutations_table.py",
        get_script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/sql/get_mutations.sql",
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"]
    conda:
        CFG["conda_envs"]["hotmaps"]
    resources:
        **CFG["resources"]["get_mutations"]
    shell:
        op.as_one_line("""
        python {params.load_script} -m {input.mysql_mut} 
        --host {params.mysql_host} --mysql-user {params.mysql_user} 
        --mysql-passwd {params.mysql_pass} --db {params.mysql_db} 
        && 
        mysql -u {params.mysql_user} -A -p{params.mysql_pass} 
        -h {params.mysql_host} {params.mysql_db} < {params.get_script} > {output.mutations}
        """)

rule _hotmaps_split_pdbs:
    input:
        mutations = str(rules._hotmaps_get_mutations.output.mutations),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)
    output:
        done = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/split_pdbs.success"
    params:
        splits = CFG["options"]["hotmaps"]["pdb_splits"],
        split_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/divide_pdb_info.py"
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/split_pdbs.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/split_pdbs.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        -f {input.pdb} 
        -m {input.mutations} 
        -n {params.splits} 
        --split-dir {params.split_dir} 
        > {log.stdout} 2> {log.stderr} &&
        touch {output.done}
        """)

rule _hotmaps_run_hotspot:
    input:
        mutations = str(rules._hotmaps_split_pdbs.output.done)
    output:
        hotspot = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot/full_output/output_{split}.txt"
    params:
        script =  CFG["dirs"]["inputs"] + "HotMAPS-master/hotspot.py",
        mutation = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/mut_info_split_{split}.txt",
        pdb = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/split_pdbs/pdb_info_split_{split}.txt",
        ttype = lambda w: w.sample_set,
        num_sims = CFG["options"]["hotmaps"]["num_sims"],
        radius = CFG["options"]["hotmaps"]["radius"],
        stop_criteria = CFG["options"]["hotmaps"]["stop_criteria"],
        error = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot/error/error_pdb_{split}.txt"
    conda: 
        CFG["conda_envs"]["hotmaps"]
    threads: CFG["threads"]["hotmaps"]
    resources:
        **CFG["resources"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/run_hotspot/hotspot_{split}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/run_hotspot/hotspot_{split}.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        --log-level=INFO 
        -m {params.mutation} 
        -a {params.pdb} 
        -t {params.ttype} 
        -n {params.num_sims} 
        -r {params.radius} 
        -o {output.hotspot} 
        -e {params.error} 
        --log {log.stdout} 
        2> {log.stderr}
        """)

def _get_splits(wildcards):
    CFG = config["lcr-modules"]["hotmaps"]
    SPLITS = []
    for i in range(0, CFG["options"]["hotmaps"]["pdb_splits"]):
        SPLITS.append(str(i))
    return expand(str(rules._hotmaps_run_hotspot.output.hotspot), split = SPLITS, allow_missing=True)

rule _hotmaps_merge_hotspots:
    input:
        mutations = _get_splits
    output:
        merged = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/output_merged.txt"
    params:
        output_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot/full_output/output_*"
    shell:
        op.as_one_line("""
        cat {params.output_dir} | 
        awk -F"\t" '/Structure/{{s++}}{{if(s==1){{print$0}}}}{{if(s>1 && $0 !~ "^Structure"){{print $0}}}}' 
        > {output.merged}
        """)

rule _hotmaps_multiple_test_correct:
    input:
        merged = str(rules._hotmaps_merge_hotspots.output.merged),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation)
    output:
        mtc = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mtc_output_{q_value}.txt",
        significance = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/significance_level_{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/multiple_testing_correction.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        group_func = CFG["options"]["hotmaps"]["group_func"],
        q_value = lambda w: w.q_value,
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/multiple_test_correct/mtc_output_{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/multiple_test_correct/mtc_output_{q_value}.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        -i {input.merged} 
        -f {params.group_func} 
        -m {params.mupit_dir} 
        -q {params.q_value} 
        -o {output.mtc} 
        -s {output.significance} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_find_gene:
    input:
        mtc = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation)
    output:
        hotspots = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot_regions_gene_{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/find_hotspot_regions_gene.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/find_gene/find_gene_{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/find_gene/find_gene_{q_value}.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        -m {input.mtc} 
        -a {params.mupit_dir} 
        -p {input.pdb} 
        -r {params.radius} 
        -q {params.q_value} 
        -o {output.hotspots} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_find_structure:
    input:
        merged = str(rules._hotmaps_merge_hotspots.output.merged),
        mtc = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation),
        significance = str(rules._hotmaps_multiple_test_correct.output.significance)
    output:
        structures = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot_regions_structure_{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/find_hotspot_regions_struct.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    conda:
        CFG["conda_envs"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/find_structure/find_structure_{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/find_structure/find_structure_{q_value}.stderr.log"
    shell:
        op.as_one_line("""
        python {params.script} 
        -i {input.merged} 
        -a {params.mupit_dir} 
        -p {input.pdb} 
        -r {params.radius} 
        -o {output.structures} 
        -s {input.significance} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_detailed_hotspots:
    input:
        hotspots = str(rules._hotmaps_find_gene.output.hotspots),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation),
        merged_output = str(rules._hotmaps_merge_hotspots.output.merged),
        mtc_file = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb_info = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)
    output:
        coordinates = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/genomic_coordinates_hotspot_regions_gene_{q_value}.txt",
        detailed = CFG["dirs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/detailed_hotspot_regions_gene_{q_value}.txt"
    params:
        script = CFG["detailed_hotspots_script"],
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    log:
        stdout = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/detailed_hotspots/detailed_hotspots_{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/detailed_hotspots/detailed_hotspots_{q_value}.stderr.log"
    threads: CFG["threads"]["hotmaps"]
    resources:
        **CFG["resources"]["hotmaps"]
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        {params.script} 
        --hotspots {input.hotspots} 
        --mupit-annotation {input.mupit_annotation} 
        --output-merged {input.merged_output} 
        --mtc-file {input.mtc_file} 
        --q-value {params.q_value} 
        --pdb-info {input.pdb_info} 
        --angstroms {params.radius} 
        --coordinates-out {output.coordinates} 
        --enriched-out {output.detailed} 
        --overwrite > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_output:
    input:
        hotspots = str(rules._hotmaps_find_gene.output.hotspots),
        structures = str(rules._hotmaps_find_structure.output.structures),
        detailed = str(rules._hotmaps_detailed_hotspots.output.detailed),
        coordinates = str(rules._hotmaps_detailed_hotspots.output.coordinates)
    output:
        hotspots = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot_regions_gene_{q_value}.txt",
        structures = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/hotspot_regions_struct_{q_value}.txt",
        detailed = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/detailed_hotspot_regions_gene_{q_value}.txt",
        coordinates = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/{md5sum}/genomic_coordinates_hotspot_regions_gene_{q_value}.txt"
    run:
        op.relative_symlink(input.hotspots, output.hotspots, in_module=True)
        op.relative_symlink(input.structures, output.structures, in_module=True)
        op.relative_symlink(input.detailed, output.detailed, in_module = True)
        op.relative_symlink(input.coordinates, output.coordinates, in_module = True)

def _for_aggregate(wildcards):
    CFG = config["lcr-modules"]["hotmaps"]
    checkpoint_output = os.path.dirname(str(checkpoints._hotmaps_prep_input.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output+"/{md5sum}.maf.content")
    return expand(
        [
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/hotspot_regions_gene_{{q_value}}.txt",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/hotspot_regions_struct_{{q_value}}.txt",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/detailed_hotspot_regions_gene_{{q_value}}.txt",
            CFG["dirs"]["outputs"] + "{{genome_build}}/{{sample_set}}--{{launch_date}}/{md5sum}/genomic_coordinates_hotspot_regions_gene_{{q_value}}.txt"
        ],
        md5sum = SUMS
    )

rule _hotmaps_aggregate:
    input:
        _for_aggregate
    output:
        aggregate = CFG["dirs"]["outputs"] + "{genome_build}/{sample_set}--{launch_date}/aggregate/{sample_set}--{launch_date}--q{q_value}.done"
    shell:
        "touch {output.aggregate}"

# Generates the target sentinels for each run, which generate the symlinks
rule _hotmaps_all:
    input:
        expand(
            [
                CFG["dirs"]["inputs"] + "maf/{genome_build}/{sample_set}--{launch_date}/done",
                str(rules._hotmaps_aggregate.output.aggregate)
            ],
            genome_build = "grch37",
            sample_set = CFG["maf_processing"]["sample_sets"],
            launch_date = launch_date,
            q_value = CFG["options"]["hotmaps"]["q_value"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
