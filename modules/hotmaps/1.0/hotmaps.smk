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
        modbase_dir = "modbase_dir=" + CFG["pdb_structure_dirs"]["modbase_dir"],
        pdb_dir = "pdb_dir=" + CFG["pdb_structure_dirs"]["pdb_dir"],
        refseq_homology_dir = "refseq_homology: " + CFG["pdb_structure_dirs"]["refseq_homology_dir"],
        ensembl_homology_dir = "ensembl_homology: " + CFG["pdb_structure_dirs"]["ensembl_homology_dir"],
        biological_assembly_dir = "biological_assembly: " + CFG["pdb_structure_dirs"]["biological_assembly_dir"],
        non_biological_assembly_dir = "non_biological_assembly: " + CFG["pdb_structure_dirs"]["non_biological_assembly_dir"],
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
rule _hotmaps_input_maf:
    input:
        maf = CFG["inputs"]["input_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}/input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)

rule _hotmaps_input_subsets:
    input:
        sample_subsets = CFG["inputs"]["sample_sets"]
    output:
        sample_subsets = CFG["dirs"]["inputs"] + "sample_sets/sample_sets.tsv"
    run:
        op.absolute_symlink(input.sample_subsets, output.sample_subsets)

rule _hotmaps_prep_input:
    input:
        maf = expand(
            str(rules._hotmaps_input_maf.output.maf),
            allow_missing=True,
            seq_type = CFG["maf_processing"]["seq_types"]
            ),
        sample_sets = ancient(str(rules._hotmaps_input_subsets.output.sample_subsets))
    output:
        #maf_hotmaps = CFG["dirs"]["inputs"] + "maf/input.{sample_set}.maf",
        maf = CFG["dirs"]["inputs"] + "maf/{sample_set}/{sample_set}.maf"
    log:
        stdout = CFG["logs"]["inputs"] + "{sample_set}/prepare_maf.stdout.log",
        stderr = CFG["logs"]["inputs"] + "{sample_set}/prepare_maf.stderr.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["maf_processing"]["include_non_coding"]).upper(),
        script = CFG["maf_processing"]["prepare_mafs"]
    shell:
        op.as_one_line("""
        Rscript {params.script}
        {input.maf}
        {input.sample_sets}
        $(dirname {output.maf})/
        {wildcards.sample_set}
        HotMAPS
        {params.include_non_coding}
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_split_dnps:
    input:
        maf = str(rules._hotmaps_prep_input.output.maf)
    output:
        dnps = CFG["dirs"]["inputs"] + "maf2vcf/{sample_set}/maf/{sample_set}.dnps.maf",
        filtered_maf = CFG["dirs"]["inputs"] + "maf2vcf/{sample_set}/maf/{sample_set}.dnp_filtered.maf"
    shell:
        op.as_one_line("""
        variant_type_col=$(head -n 1 {input.maf} | sed 's/\\t/\\n/g' | nl | grep "Variant_Type" | cut -f 1) &&
        protein_position_col=$(head -n 1 {input.maf} | sed 's/\\t/\\n/g' | nl | grep "Protein_position" | cut -f 1) &&
        cat <( head -n 1 {input.maf} ) <( awk -v var_col="$variant_type_col" -v protein_col="$protein_position_col" ' {{ if ( $var_col=="DNP" && $protein_col ~ /[0-9?]+-[0-9?]+/) print $0 }} ' {input.maf} ) > {output.dnps} &&
        awk -v var_col="$variant_type_col" -v protein_col="$protein_position_col" ' {{ if ( $var_col != "DNP" || $protein_col !~ /[0-9?]+-[0-9?]+/ ) print $0 }} ' {input.maf} > {output.filtered_maf}
        """)

checkpoint _hotmaps_maf2vcf:
    input:
        dnps = str(rules._hotmaps_split_dnps.output.dnps),
        fasta = reference_files("genomes/grch37/genome_fasta/genome.fa")
    output:
        vcf = CFG["dirs"]["inputs"] + "maf2vcf/{sample_set}/vcf/{sample_set}.dnps.vcf"
    params:
        vcf_dir = CFG["dirs"]["inputs"] + "maf2vcf/{sample_set}/vcf"
    conda:
        CFG["conda_envs"]["bcftools"]
    log:
        stdout = CFG["logs"]["inputs"] + "maf2vcf/{sample_set}/{sample_set}.maf2vcf.stdout.log",
        stderr = CFG["logs"]["inputs"] + "maf2vcf/{sample_set}/{sample_set}.maf2vcf.stderr.log"
    shell:
        op.as_one_line("""
        maf2vcf.pl 
        --input-maf {input.dnps} 
        --output-dir {params.vcf_dir} 
        --output-vcf {output.vcf} 
        --ref-fasta {input.fasta} 
        --per-tn-vcfs 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_bcftools:
    input:
        vcf = CFG["dirs"]["inputs"] + "maf2vcf/{sample_set}/vcf/{tumour_id}_vs_{normal_sample_id}.vcf",
        dnp_vcf = str(rules._hotmaps_maf2vcf.output.vcf)
    output:
        vcf = CFG["dirs"]["inputs"] + "bcftools/{sample_set}/vcf/{tumour_id}_vs_{normal_sample_id}.annotate.vcf"
    conda:
        CFG["conda_envs"]["bcftools"] 
    log:
        stdout = CFG["logs"]["inputs"] + "bcftools/{sample_set}/{tumour_id}_vs_{normal_sample_id}.annotate.log",
        stderr = CFG["logs"]["inputs"] + "bcftools/{sample_set}/{tumour_id}_vs_{normal_sample_id}.annotate.log"
    shell:
        op.as_one_line("""
        bcftools norm --atomize 
        --output {output.vcf} 
        --output-type v {input.vcf} 
        > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_vcf2maf:
    input:
        vcf = str(rules._hotmaps_bcftools.output.vcf),
        fasta = reference_files("genomes/grch37/genome_fasta/genome.fa"),
        vep_cache = CFG["vcf2maf"]["vep_cache"]
    output:
        maf = CFG["dirs"]["inputs"] + "vcf2maf/{sample_set}/maf/{tumour_id}_vs_{normal_sample_id}.maf",
        vep = CFG["dirs"]["inputs"] + "bcftools/{sample_set}/vcf/{tumour_id}_vs_{normal_sample_id}.annotate.vep.vcf"
    params:
        opts = CFG["vcf2maf"]["options"]
    conda:
        CFG["conda_envs"]["bcftools"]
    log:
        stdout = CFG["logs"]["inputs"] + "vcf2maf/{sample_set}/{tumour_id}_vs_{normal_sample_id}.stdout.log",
        stderr = CFG["logs"]["inputs"] + "vcf2maf/{sample_set}/{tumour_id}_vs_{normal_sample_id}.stderr.log"
    shell:
        op.as_one_line("""
        vepPATH=$(dirname $(which variant_effect_predictor.pl))/../share/variant-effect-predictor* ;
        vcf2maf.pl 
        --input-vcf {input.vcf} 
        --output-maf {output.maf} 
        --tumor-id {wildcards.tumour_id} 
        --normal-id {wildcards.normal_sample_id} 
        --vcf-tumor-id TUMOR 
        --vcf-normal-id NORMAL 
        --vep-path $vepPATH 
        --vep-data {input.vep_cache}
        --ref-fasta {input.fasta}
        --retain-info gnomADg_AF,gnomADg_AF
        > {log.stdout} 2> {log.stderr}
        """) 

def _get_annotated_mafs(wildcards):
    CFG = config["lcr-modules"]["hotmaps"]
    checkpoint_outputs = checkpoints._hotmaps_maf2vcf.get(**wildcards).output.vcf

    # Get rows of DNPs that were reannotated
    dnps = expand(
        str(rules._hotmaps_split_dnps.output.dnps),
        sample_set = wildcards.sample_set
    )

    # Read in MAF as table
    if os.path.exists(dnps[0]):
        dnp_table = pd.read_table(dnps[0], comment="#", sep="\t")
        sample_table = dnp_table[["Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode"]].drop_duplicates()
        return expand(
            expand(
                str(rules._hotmaps_vcf2maf.output.maf),
                zip,
                tumour_id = sample_table["Tumor_Sample_Barcode"],
                normal_sample_id = sample_table["Matched_Norm_Sample_Barcode"],
                allow_missing=True
            ),
            sample_set = wildcards.sample_set
        )
    else:
        return []

rule _hotmaps_merge_mafs:
    input:
        maf_annotated = _get_annotated_mafs,
        maf = str(rules._hotmaps_split_dnps.output.filtered_maf)
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{sample_set}/input.{sample_set}.maf"
    run:
        print(input.maf_annotated)
        #annotated = (input.maf_annotated).split(" ")
        main_maf = pd.read_table(input.maf, comment = "#", sep="\t")
        for maf in input.maf_annotated:
            print(maf)
            df = pd.read_table(maf, comment="#", sep="\t")
            main_maf = pd.concat([main_maf, df])
        main_maf = main_maf[["Tumor_Sample_Barcode","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2","Hugo_Symbol","Variant_Classification","HGVSp_Short","Transcript_ID","Strand"]]
        main_maf.to_csv(output.maf, sep="\t", index=False)

rule _hotmaps_prep_mutations:
    input:
        maf = CFG["inputs"]["custom_maf"],
    output:
        mupit_non_filtered = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/non_filtered_mupit.input.{sample_set}.maf"
    params:
        mut_dir = lambda w: config["lcr-modules"]["hotmaps"]["dirs"]["hotmaps"] + f"/{w.sample_set}/mutations/",
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"],
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mut_regex = "'^input.+\.maf'$",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/map_maf_to_structure.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} --data-dir {params.mut_dir} --match-regex {params.mut_regex} --host {params.mysql_host} --db {params.mysql_db} --mysql-user {params.mysql_user} --mysql-passwd {params.mysql_pass} --output-dir {params.mut_dir}
        """)

rule _hotmaps_prep_mupit_annotation:
    input:
        maf = CFG["inputs"]["custom_maf"]
    output:
        annotation = CFG["dirs"]["hotmaps"] + "{sample_set}/mupit_annotations/mupit_mutations_{sample_set}"
    params:
        mut_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/",
        tumor_type = "{sample_set}",
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"],
        cov_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/",
        hypermut = "-mt " + CFG["options"]["hotmaps"]["hypermut_threshold"] if CFG["options"]["hotmaps"]["hypermut_threshold"] is not None else "",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/maf/convert_maf_to_mupit.py"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} --maf {input.maf} -mh {params.mysql_host} -mdb {params.mysql_db} --mysql-user {params.mysql_user} --mysql-passwd {params.mysql_pass} --tumor-type {params.tumor_type} --no-stratify {params.hypermut} -i {params.cov_dir} --output {output.annotation}
        """)

rule _hotmaps_filter_hypermutated:
    input:
        maf = CFG["inputs"]["custom_maf"],
        nf_mupit = str(rules._hotmaps_prep_mutations.output.mupit_non_filtered)
    output:
        hypermut = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/hypermutated.{sample_set}.txt",
        mupit = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/mupit.input.{sample_set}.maf"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/filter_hypermutated.py",
        mut_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/",
        mut_regex = "'^input.+\.maf'$",
        mut_threshold = "--mut-threshold " + CFG["options"]["hotmaps"]["hypermut_threshold"] if CFG["options"]["hotmaps"]["hypermut_threshold"] is not None else "",
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} --raw-dir {params.mut_dir} --match-regex {params.mut_regex} {params.mut_threshold} --sample-col Tumor_Sample_Barcode --data-dir {params.mut_dir}
        """)


rule _hotmaps_count_mutations:
    input:
        mupit = str(rules._hotmaps_filter_hypermutated.output.mupit)
    output:
        collected = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/collected.input.{sample_set}.maf"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/count_mutations.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/"
    shell:
        op.as_one_line("""
        python {params.script} --data-dir {params.data_dir}
        """)

rule _hotmaps_format_mutations:
    input:
        collected = str(rules._hotmaps_count_mutations.output.collected)
    output:
        tcga = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/mutation_tcga.{sample_set}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/format_mutations_table.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/"
    shell:
        op.as_one_line("""
        python {params.script} --data-dir {params.data_dir}
        """)

rule _hotmaps_merge_mutations:
    input:
        tcga = str(rules._hotmaps_format_mutations.output.tcga)
    output:
        mysql_mut = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/mysql.mutations.tcga.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/merge_mutations_table_data.py",
        data_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/"
    shell:
        op.as_one_line("""
        python {params.script} {params.data_dir}
        """)

rule _hotmaps_load_mutations:
    input:
        mysql_mut = str(rules._hotmaps_merge_mutations.output.mysql_mut)
    output:
        mysql_loaded = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/loaded.success"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/mupit/load_mutations_table.py",
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"]
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} -m {input.mysql_mut} --host {params.mysql_host} --mysql-user {params.mysql_user} --mysql-passwd {params.mysql_pass} --db {params.mysql_db} &&
        touch {output.mysql_loaded}
        """)

rule _hotmaps_get_mutations:
    input:
        done = str(rules._hotmaps_load_mutations.output.mysql_loaded)
    output:
        mutations = CFG["dirs"]["hotmaps"] + "{sample_set}/mutations/mutations.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/sql/get_mutations.sql",
        mysql_db = CFG["options"]["mysql"]["mysql_db"],
        mysql_host = CFG["options"]["mysql"]["mysql_host"],
        mysql_user = CFG["options"]["mysql"]["mysql_user"],
        mysql_pass = CFG["options"]["mysql"]["mysql_passwd"]
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        mysql -u {params.mysql_user} -A -p{params.mysql_pass} -h {params.mysql_host} {params.mysql_db} < {params.script} > {output.mutations}
        """)

rule _hotmaps_split_pdbs:
    input:
        mutations = str(rules._hotmaps_get_mutations.output.mutations),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)
    output:
        done = CFG["dirs"]["hotmaps"] + "{sample_set}/split_pdbs/split_pdbs.success"
    params:
        splits = CFG["options"]["hotmaps"]["pdb_splits"],
        split_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/split_pdbs/",
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/scripts/divide_pdb_info.py"
    shell:
        op.as_one_line("""
        python {params.script} -f {input.pdb} -m {input.mutations} -n {params.splits} --split-dir {params.split_dir} &&
        touch {output.done}
        """)

rule _hotmaps_run_hotspot:
    input:
        mutations = str(rules._hotmaps_split_pdbs.output.done)
    output:
        hotspot = CFG["dirs"]["hotmaps"] + "{sample_set}/hotspot/full_output/output_{split}.txt"
    params:
        script =  CFG["dirs"]["inputs"] + "HotMAPS-master/hotspot.py",
        mutation = CFG["dirs"]["hotmaps"] + "{sample_set}/split_pdbs/mut_info_split_{split}.txt",
        pdb = CFG["dirs"]["hotmaps"] + "{sample_set}/split_pdbs/pdb_info_split_{split}.txt",
        num_sims = CFG["options"]["hotmaps"]["num_sims"],
        radius = CFG["options"]["hotmaps"]["radius"],
        error = CFG["dirs"]["hotmaps"] + "{sample_set}/hotspot/error/error_pdb_{split}.txt"
    conda: 
        CFG["conda_envs"]["hotmaps"]
    threads: CFG["threads"]["hotmaps"]
    resources:
        **CFG["resources"]["hotmaps"]
    log:
        stdout = CFG["logs"]["hotmaps"] + "{sample_set}/hotmaps_run_{split}.stdout.log"
    shell:
        op.as_one_line("""
        python {params.script} --log-level=INFO -m {params.mutation} -a {params.pdb} -t EVERY -n {params.num_sims} 
        -r {params.radius} -o {output.hotspot} -e {params.error} --log {log.stdout}
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
        merged = CFG["dirs"]["hotmaps"] + "{sample_set}/output_merged.txt"
    params:
        output_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/hotspot/full_output/output_*"
    shell:
        op.as_one_line("""
        cat {params.output_dir} | awk -F"\t" '/Structure/{{s++}}{{if(s==1){{print$0}}}}{{if(s>1 && $0 !~ "^Structure"){{print $0}}}}' > {output.merged}
        """)

rule _hotmaps_multiple_test_correct:
    input:
        merged = str(rules._hotmaps_merge_hotspots.output.merged),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation)
    output:
        mtc = CFG["dirs"]["hotmaps"] + "{sample_set}/mtc_output_min_.{q_value}.txt",
        significance = CFG["dirs"]["hotmaps"] + "{sample_set}/significance_level_.{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/multiple_testing_correction.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        group_func = CFG["options"]["hotmaps"]["group_func"],
        q_value = lambda w: w.q_value,
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} -i {input.merged} -f {params.group_func} -m {params.mupit_dir} -q {params.q_value} -o {output.mtc} -s {output.significance}
        """)

rule _hotmaps_find_gene:
    input:
        mtc = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation)
    output:
        hotspots = CFG["dirs"]["hotmaps"] + "{sample_set}/hotspot_regions_gene_.{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/find_hotspot_regions_gene.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} -m {input.mtc} -a {params.mupit_dir} -p {input.pdb} -r {params.radius} -q {params.q_value} -o {output.hotspots}
        """)

rule _hotmaps_find_structure:
    input:
        merged = str(rules._hotmaps_merge_hotspots.output.merged),
        mtc = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation),
        significance = str(rules._hotmaps_multiple_test_correct.output.significance)
    output:
        structures = CFG["dirs"]["hotmaps"] + "{sample_set}/hotspot_regions_structure_.{q_value}.txt"
    params:
        script = CFG["dirs"]["inputs"] + "HotMAPS-master/find_hotspot_regions_struct.py",
        mupit_dir = CFG["dirs"]["hotmaps"] + "{sample_set}/mupit_annotations/",
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        python {params.script} -i {input.merged} -a {params.mupit_dir} -p {input.pdb} -r {params.radius} -o {output.structures} -s {input.significance}
        """)

rule _hotmaps_get_detailed_hotspots:
    input:
        hotspots = str(rules._hotmaps_find_gene.output.hotspots),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation),
        merged_output = str(rules._hotmaps_merge_hotspots.output.merged),
        mtc_file = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb_info = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)
    output:
        hotspots = CFG["dirs"]["hotmaps"] + "{sample_set}/enriched_hotspot_regions_gene_.{q_value}.txt"
    params:
        script = CFG["options"]["enrich_hotmaps_script"],
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    log:
        stdout = CFG["logs"]["hotmaps"] + "{sample_set}/get_detailed_hotspots_.{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{sample_set}/get_detailed_hotspots.{q_value}.stderr.log"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        {params.script} --hotspots {input.hotspots} --mupit-annotation {input.mupit_annotation} 
        --output-merged {input.merged_output} --mtc-file {input.mtc_file} 
        --q-value {params.q_value} --pdb-info {input.pdb_info} 
        --angstroms {params.radius} --metadata-out {output.hotspots} > {log.stdout} 2> {log.stderr}
        """)

rule _hotmaps_get_hotspot_coordinates:
    input:
        hotspots = str(rules._hotmaps_find_gene.output.hotspots),
        mupit_annotation = str(rules._hotmaps_prep_mupit_annotation.output.annotation),
        merged_output = str(rules._hotmaps_merge_hotspots.output.merged),
        mtc_file = str(rules._hotmaps_multiple_test_correct.output.mtc),
        pdb_info = str(rules._hotmaps_add_pdb_description.output.fully_described_pdb)
    output:
        hotspots = CFG["dirs"]["hotmaps"] + "{sample_set}/genomic_coordinates_hotspot_regions_gene_.{q_value}.txt"
    params:
        script = CFG["options"]["enrich_hotmaps_script"],
        radius = CFG["options"]["hotmaps"]["radius"],
        q_value = lambda w: w.q_value
    log:
        stdout = CFG["logs"]["hotmaps"] + "{sample_set}/get_hotspot_coordinates.{q_value}.stdout.log",
        stderr = CFG["logs"]["hotmaps"] + "{sample_set}/get_hotspot_coordinates.{q_value}.stderr.log"
    conda:
        CFG["conda_envs"]["hotmaps"]
    shell:
        op.as_one_line("""
        {params.script} --hotspots {input.hotspots} --mupit-annotation {input.mupit_annotation} 
        --output-merged {input.merged_output} --mtc-file {input.mtc_file} 
        --q-value {params.q_value} --pdb-info {input.pdb_info} 
        --angstroms {params.radius} --metadata-out {output.hotspots} --maf-mode > {log.stdout} 2> {log.stderr}
        """)


rule _hotmaps_symlink_output:
    input:
        hotspots = str(rules._hotmaps_find_gene.output.hotspots),
        structures = str(rules._hotmaps_find_structure.output.structures),
        enriched = str(rules._hotmaps_get_detailed_hotspots.output.hotspots),
        coordinates = str(rules._hotmaps_get_hotspot_coordinates.output.hotspots)
    output:
        hotspots = CFG["dirs"]["outputs"] + "{sample_set}/hotspot_regions_gene_.{q_value}.txt",
        structures = CFG["dirs"]["outputs"] + "{sample_set}/hotspot_regions_struct_.{q_value}.txt",
        enriched = CFG["dirs"]["outputs"] + "{sample_set}/enriched_hotspot_regions_gene_.{q_value}.txt",
        coordinates = CFG["dirs"]["outputs"] + "{sample_set}/genomic_coordinates_hotspot_regions_gene_.{q_value}.txt"
    run:
        op.relative_symlink(input.hotspots, output.hotspots, in_module=True)
        op.relative_symlink(input.structures, output.structures, in_module=True)
        op.relative_symlink(input.enriched, output.enriched, in_module = True)
        op.relative_symlink(input.coordinates, output.coordinates, in_module = True)


# Generates the target sentinels for each run, which generate the symlinks
rule _hotmaps_all:
    input:
        expand(
            [
                str(rules._hotmaps_symlink_output.output.hotspots),
                str(rules._hotmaps_symlink_output.output.structures),
                str(rules._hotmaps_symlink_output.output.enriched),
                str(rules._hotmaps_symlink_output.output.coordinates)
            ],
            sample_set = CFG["maf_processing"]["sample_sets"],
            q_value = CFG["options"]["hotmaps"]["q_value"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
