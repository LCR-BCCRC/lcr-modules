#!/usr/bin/env snakemake
import os
import pandas
import glob as glob

##### ATTRIBUTION #####


# Original Author:  Chris Rushton
# Module Author:    Chris Rushton
# Contributors:     NA

# LymphGen classifier writen by George Wright
# DO NOT DISTRIBUTE THE LYMPHGEN SOURCE CODE WITHOUT GEORGE WRIGHT'S PERMISSION
# This module has been brought to you by Krysta's Girl Guide cookies. Helping me
# miss all weight loss goals since 2018

# This module can only be run with grch37/hg19 genome builds

##### SETUP #####


# Import package with useful functions for developing analysis modules
import oncopipe as op
#
assert (config["lcr-modules"]["lymphgen"]["inputs"]["sample_seg"]["empty"] == "{MODSDIR}/etc/empty.seg"), (
    "'config['inputs']['sample_seg']['empty']' must be present and unmodified"
)

assert (config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["other"]["empty_sv"] == "{MODSDIR}/etc/empty_svar_master.tsv"), (
    "'config['inputs']['sample_sv_info']['other']['empty_sv']' must be present and unmodified"
)

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["lymphgen"]`
CFG = op.setup_module(
    name = "lymphgen",
    version = "2.1",
    subdirectories = ["inputs", "reformat_seg", "lymphgen_input", "lymphgen_run", "composite_other", "outputs"],
)


def _find_best_seg(wildcards):
    this_tumor = op.filter_samples(RUNS, tumour_genome_build = wildcards.genome_build, tumour_sample_id = wildcards.tumour_id, normal_sample_id = wildcards.normal_id, pair_status = wildcards.pair_status, tumour_seq_type = wildcards.seq_type)
    return this_tumor.cnv_path

def _find_best_sv(wildcards):
    this_tumor = op.filter_samples(RUNS, tumour_genome_build = wildcards.genome_build, tumour_sample_id = wildcards.tumour_id, normal_sample_id = wildcards.normal_id, pair_status = wildcards.pair_status, tumour_seq_type = wildcards.seq_type)
    if this_tumor["has_sv"].bool() == False:
        return config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["other"]["empty_sv"]
    else:
        return this_tumor.sv_path

def has_fish(mod_config):
    fish = CFG["inputs"]["sample_sv_info"]["fish"]
    fish = pandas.read_csv(fish, sep='\t')
    fish_sample_ids = fish['sample_id'].to_numpy()
    mod_config['has_fish'] = False
    for i, row in mod_config.iterrows():
        if row['tumour_sample_id'] in fish_sample_ids:
            mod_config['has_fish'] = True
        else:
            mod_config['has_fish'] = False

def best_seg_sv(mod_config):
    # set the cnv_path to be used for each sample based on sample_seg keys in the config 
    callers_seg = list(mod_config["inputs"]["sample_seg"].keys())
    callers_seg = [caller.lower() for caller in callers_seg]
    dfm = mod_config["runs"]
    df = dfm.to_dict('index')
    for index in df.keys():
        possible_outputs = []
        for caller in callers_seg:
            if caller == "empty":
                possible_output = [mod_config["inputs"]["sample_seg"][caller]]
            else:
                possible_output = (expand(mod_config["inputs"]["sample_seg"][caller], zip,
                                        tumour_id=df[index]["tumour_sample_id"],
                                        normal_id=df[index]["normal_sample_id"],
                                        pair_status=df[index]["pair_status"],
                                        genome_build=df[index]["tumour_genome_build"],
                                        seq_type=df[index]["tumour_seq_type"]))
        
            if os.path.exists(possible_output[0]):
                possible_outputs += possible_output
            
        df[index]["cnv_path"] = possible_outputs[0]
        
        mod_config["runs"].at[index, "cnv_path"] = df[index]["cnv_path"]

    # set the sv_path to be used for each sample based on ["sample_sv_info"]["other"] keys 
    callers_sv = list(mod_config["inputs"]["sample_sv_info"]["other"].keys())
    callers_sv = [caller.lower() for caller in callers_sv]
    dfm = mod_config["runs"]
    df = dfm.to_dict('index')
    for index in df.keys():
        possible_outputs = []
        for caller in callers_sv:
            if caller == "empty_sv":
                possible_output = [mod_config["inputs"]["sample_sv_info"]["other"][caller]]
            else:
                possible_output = (expand(mod_config["inputs"]["sample_sv_info"]["other"][caller], zip,
                                        tumour_id=df[index]["tumour_sample_id"],
                                        normal_id=df[index]["normal_sample_id"],
                                        pair_status=df[index]["pair_status"],
                                        genome_build=df[index]["tumour_genome_build"],
                                        seq_type=df[index]["tumour_seq_type"]))

            if os.path.exists(possible_output[0]):
                possible_outputs += possible_output

        df[index]["sv_path"] = possible_outputs[0]

        mod_config["runs"].at[index, "sv_path"] = df[index]["sv_path"]

    all = mod_config["runs"]

    # creates a column in the RUNS table for having sv 
    # if a sample is capture, it should only go through the with_sv pipeline if it has FISH
    #       this is because detecting SVs using short read sequencing often resembles common #       sequencing and alignment artifacts. 
    # if a sample is genome, it should go through with_sv whether it has FISH OR svar_master OR # manta outputs
    all['has_sv'] = False
    for i, row in all.iterrows():
        # if the sample is genome and it doesn't have an empty SV path for manta or svar_master 
        if (row['tumour_seq_type'] == "genome") & (row['sv_path'] != config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["other"]["empty_sv"]):
            all.at[i, 'has_sv'] = True
        else:
            all.at[i, 'has_sv'] = False

    # this function creates a column 'has_fish' based on whether there is FISH data for the 
    # sample in the input FISH table 
    has_fish(all)

    return(all)


RUNS = best_seg_sv(CFG)


# Define rules to be run locally when using a compute cluster
# I put everything under this, since these rules don't take very long
localrules:
    _install_lgenic,
    _lymphgen_input_cnv,
    _lymphgen_input_no_cnv,
    _lymphgen_gamblr_config,
    _lymphgen_add_sv,
    _lymphgen_add_sv_blank,
    _lymphgen_run_cnv_A53,
    _lymphgen_run_cnv_noA53, 
    _lymphgen_run_no_cnv,
    _lymphgen_reformat_seg,
    _lymphgen_output_txt,
    _lymphgen_all,


# Sanitize input
lgenic_path = CFG["inputs"]["lgenic_exec"]
if not lgenic_path.endswith(os.sep) and lgenic_path != "":
    CFG["inputs"]["lgenic_exec"] = CFG["inputs"]["lgenic_exec"] + os.sep


##### RULES #####


# DOWNLOAD CHRIS'S VERY AWESOME LYMPHGEN CONVERSION SCRIPT. ALL CAPS
rule _install_lgenic:
    params:
        lgenic_dir = CFG["inputs"]["lgenic_exec"]
    output:
        lgenic_script = CFG["inputs"]["lgenic_exec"] + "generate_input.py",
        lymphgen_genes = CFG["inputs"]["lgenic_exec"] + "resources" + os.sep + "lymphgen_genes.txt",
        hugo2entrez = CFG["inputs"]["lgenic_exec"] + "resources" + os.sep + "hugo2entrez.tsv",
        gene_coords = CFG["inputs"]["lgenic_exec"] + "resources" + os.sep + "gene_coordinates.GRCh37.bed6",
        arm_coords = CFG["inputs"]["lgenic_exec"] + "resources" + os.sep + "chrom_arm.hg19.tsv"
    group:
        "lymphgen"
    shell:
        '''
        download_url=$(curl --silent "https://api.github.com/repos/LCR-BCCRC/LGenIC/releases/169598665" | grep 'tarball_url' | sed 's/.*:[ ]//' | sed 's/,$//' | sed 's/"//g');
        mkdir -p {params.lgenic_dir};

        wget -cO - $download_url > {params.lgenic_dir}/LGenIC.tar.gz && tar -C {params.lgenic_dir} -xf {params.lgenic_dir}/LGenIC.tar.gz && rm {params.lgenic_dir}/LGenIC.tar.gz;
        mv {params.lgenic_dir}/ckrushton-LGenIC-*/* {params.lgenic_dir}/ && rm -r {params.lgenic_dir}/ckrushton-LGenIC-*/;
        '''


# STEP 1: INPUT SYMLINKS
# Symlinks the input files into the module results directory (under '00-inputs/')

rule _lymphgen_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"] 
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    group:
        "lymphgen"
    run:
        op.relative_symlink(input.maf, output.maf)

rule _lymphgen_input_seg:
    input:
        seg = _find_best_seg
    output:
        seg = CFG["dirs"]["inputs"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.seg"
    group:
        "lymphgen"
    run:
        op.relative_symlink(input.seg, output.seg)
        
rule _lymphgen_input_sv: 
    input: 
        sv = _find_best_sv
    output: 
        sv = CFG["dirs"]["inputs"] + "sv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    group: 
        "lymphgen"
    run: 
        op.relative_symlink(input.sv, output.sv)


# STEP 2: REFORMAT SEG FILE
# Make sure the SEG columns are consistent
rule _lymphgen_reformat_seg:
    input:
        seg = str(rules._lymphgen_input_seg.output.seg)
    output:
        seg = CFG["dirs"]["reformat_seg"] + "seg/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}--reformat.seg"
    params:
        tumor_sample_barcode_name = CFG["options"]["reformat_seg"]["Tumor_Sample_Barcode"],
        chromosome_name = CFG["options"]["reformat_seg"]["chromosome"],
        start_name = CFG["options"]["reformat_seg"]["start"],
        end_name = CFG["options"]["reformat_seg"]["end"],
        cn_name = CFG["options"]["reformat_seg"]["CN"]
    group:
        "lymphgen"
    run:

        loaded_seg = pandas.read_csv(input.seg, sep="\t")
        seg_header = list(loaded_seg.columns)
        # Rename relevant columns
        new_cols = {
            "Tumor_Sample_Barcode": params.tumor_sample_barcode_name,
            "chromosome": params.chromosome_name,
            "start": params.start_name,
            "end": params.end_name,
            "CN": params.cn_name
        }

        for new_name, old_name in new_cols.items():
            try:
                loc = seg_header.index(old_name)
            except ValueError as e:
                raise AttributeError(f"Unable to locate column {old_name} in the SEG file {input.seg}") from e
            seg_header[loc] = new_name

        # Write out the renamed SEG file
        loaded_seg.columns = seg_header
        loaded_seg.to_csv(output.seg, sep="\t", header=True, index=False)


# STEP 3: REFORMAT INPUT TO RUN LYMPHGEN
# Reformats MAF/SEG SNV/CNV calls for LymphGen using my LGenIC script

def _get_gene_list(w): 
    CFG = config["lcr-modules"]["lymphgen"]
    this_tumor = op.filter_samples(RUNS, tumour_genome_build = w.genome_build, tumour_sample_id = w.tumour_id, normal_sample_id = w.normal_id, pair_status = w.pair_status, tumour_seq_type = w.seq_type)
    if "tumour_capture_space" in this_tumor.columns: 
        capture_space = this_tumor["tumour_capture_space"]
        capture_path = expand(CFG["inputs"]["gene_list"], capture_space = capture_space)[0]
    else: 
        capture_path = CFG["inputs"]["gene_list"]
    if os.path.exists(capture_path): 
        return str(capture_path)
    else: 
        return str(rules._install_lgenic.output.lymphgen_genes)


# With CNVs
rule _lymphgen_input_cnv:
    input:
        maf = str(rules._lymphgen_input_maf.output.maf),
        seg = str(rules._lymphgen_reformat_seg.output.seg),
        # Software and resource dependencies from LGenIC
        lgenic_script = str(rules._install_lgenic.output.lgenic_script),
        lymphgen_genes = _get_gene_list,
        hugo2entrez = str(rules._install_lgenic.output.hugo2entrez),
        gene_coords = str(rules._install_lgenic.output.gene_coords),
        arm_coords = str(rules._install_lgenic.output.arm_coords),
    output:
        sample_annotation = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_sample_annotation.tsv",
        mutation_flat = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_mutation_flat.tsv",
        gene_list = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_gene_list.txt",
        cnv_flat = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_cnv_flat.tsv",
        cnv_arm = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_cnv_arm.tsv"
    log:
        stdout = CFG["logs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}.LGenIC.stdout.log",
        stderr = CFG["logs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}.LGenIC.stderr.log"
    group:
        "lymphgen"
    params:
        outprefix = "{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}",
        logratio = "--log2" if CFG["options"]["lymphgen_input"]["use_log_ratio"].lower() == "true" else ""
    conda:
        CFG['conda_envs']['sorted_containers']
    wildcard_constraints:
        cnvs_wc = "with_cnvs"
    shell:
        op.as_one_line("""
        python {input.lgenic_script} --lymphgen_genes {input.lymphgen_genes} --outdir $(dirname {output.sample_annotation})
        --outprefix {params.outprefix} -v INFO --maf {input.maf} --entrez_ids {input.hugo2entrez} --cnvs {input.seg} {params.logratio} --genes {input.gene_coords} --arms {input.arm_coords}
        > {log.stdout} 2> {log.stderr}
        """)

# No CNVs
rule _lymphgen_input_no_cnv:
    input:
        maf = str(rules._lymphgen_input_maf.output.maf),
        # Software and resource dependencies from LGenIC
        lgenic_script = str(rules._install_lgenic.output.lgenic_script),
        lymphgen_genes = _get_gene_list,
        hugo2entrez = str(rules._install_lgenic.output.hugo2entrez),
    output:
        sample_annotation = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_sample_annotation.tsv",
        mutation_flat = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_mutation_flat.tsv",
        gene_list = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}_gene_list.txt",
    log:
        stdout = CFG["logs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}.LGenIC.stdout.log",
        stderr = CFG["logs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}.LGenIC.stderr.log"
    group:
        "lymphgen"
    params:
        outprefix = "{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}"
    conda:
        CFG['conda_envs']['sorted_containers']
    wildcard_constraints:
        cnvs_wc = "no_cnvs"
    shell:
        op.as_one_line("""
        python {input.lgenic_script} --lymphgen_genes {input.lymphgen_genes} --outdir $(dirname {output.sample_annotation})
        --outprefix {params.outprefix} -v INFO --maf {input.maf} --entrez_ids {input.hugo2entrez} > {log.stdout} 2> {log.stderr}
        """)

# STEP 4: Add SV information (if availible)

rule _lymphgen_gamblr_config:
    params:
        config_url = CFG["inputs"]["gamblr_config_url"], 
    output:
        config = "config.yml"
    group:
        "lymphgen"
    shell:
        op.as_one_line("""
        wget -qO {output.config} {params.config_url} 
        """)


rule _lymphgen_process_sv:
    input:
        fish = ancient(CFG["inputs"]["sample_sv_info"]["fish"]),
        sv = str(rules._lymphgen_input_sv.output.sv),
        gamblr = ancient(rules._lymphgen_gamblr_config.output.config)
    output:
        sv = CFG["dirs"]["inputs"] + "sv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.tsv"
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['gamblr']
    wildcard_constraints:
        sv_wc = "with_sv"
    script:
        config["lcr-modules"]["lymphgen"]['options']['lymphgen_run']['lymphgen_sv']

rule _lymphgen_add_sv:
    input:
        sample_annotation = str(rules._lymphgen_input_cnv.output.sample_annotation),
        bcl2_bcl6_sv = str(rules._lymphgen_process_sv.output.sv)
    output:
        sample_annotation = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sample_annotation.{cnvs_wc}.{sv_wc}.tsv"
    params:
        sampleIDcolname = CFG["options"]["add_svs"]["samplecol"],
        bcl2colname = CFG["options"]["add_svs"]["bcl2col"],
        bcl2true = CFG["options"]["add_svs"]["bcl2truevalues"],
        bcl2false = CFG["options"]["add_svs"]["bcl2falsevalues"],
        bcl6colname = CFG["options"]["add_svs"]["bcl6col"],
        bcl6true = CFG["options"]["add_svs"]["bcl6truevalues"],
        bcl6false = CFG["options"]["add_svs"]["bcl6falsevalues"]
    group:
        "lymphgen"
    wildcard_constraints:
        sv_wc = "with_sv"
    run:

        # Sanitize input fields
        # Since the input could be either a string or an iterable, and because python will split a string into individual letters,
        # lets make sure all input options are a set
        params.bcl2true = set([params.bcl2true] if isinstance(params.bcl2true, str) else params.bcl2true)
        params.bcl2false = set([params.bcl2false] if isinstance(params.bcl2false, str) else params.bcl2false)
        params.bcl6true = set([params.bcl6true] if isinstance(params.bcl6true, str) else params.bcl6true)
        params.bcl6false = set([params.bcl6false] if isinstance(params.bcl6false, str) else params.bcl6false)

        # Open the SV info file, and load the required columns
        loaded_sv = pandas.read_csv(input.bcl2_bcl6_sv, sep="\t")  # Pandas claims to auto-detect the file seperator, but in my experience it doesn't work
        # Check that the required columns exist
        if not params.bcl2colname in loaded_sv.columns:
            raise AttributeError("Unable to locate column \'%s\' in \'%s\'" % (params.bcl2colname, input.bcl2_bcl6_sv))
        elif not params.bcl6colname in loaded_sv.columns:
            raise AttributeError("Unable to locate column \'%s\' in \'%s\'" % (params.bcl6colname,  input.bcl2_bcl6_sv))
        elif not params.sampleIDcolname in loaded_sv.columns:
            raise AttributeError("Unable to locate column \'%s\' in \'%s\'" % (params.sampleIDcolname, input.bcl2_bcl6_sv))
        
        # Get the current sample ID
        sampleID = wildcards.tumour_id
        # Find the matching SampleID in the SV annotation file
        sampleEntry = loaded_sv.loc[loaded_sv[params.sampleIDcolname] == sampleID]
        
        # Add SV info to sample annotation file
        with open(input.sample_annotation) as f, open(output.sample_annotation, "w") as o:
            content = f.readlines()
            # write header
            header = content[0].rstrip("\n").rstrip("\r")  # Handle line endings
            o.write(header)
            o.write(os.linesep)
            
            if len(content) > 1: 
                cols = content[1].split("\t")
                try:
                    copynum = cols[1]
                except IndexError as e:
                    raise AttributeError("Input sample annotation file \'%s\' appears to be malformed" % input.sample_annotation) from e
            else: 
                copynum = "0"
                
            # BCL2
            try:
                bcl2status = sampleEntry[params.bcl2colname].tolist()[0]
                # Check to see if BCL2 is translocated or not
                if bcl2status in params.bcl2true:
                    bcl2trans = "1"
                elif bcl2status in params.bcl2false:
                    bcl2trans = "0"
                else:
                    bcl2trans = "NA"
            except IndexError:  # i.e. This sample isn't in the annotation file. set it as NA
                bcl2trans = "NA"
            # BCL6
            try:
                bcl6status = sampleEntry[params.bcl6colname].tolist()[0]
                # Check if BCL6 is translocated
                if bcl6status in params.bcl6true:
                    bcl6trans = "1"
                elif bcl6status in params.bcl6false:
                    bcl6trans = "0"
                else:
                    bcl6trans = "NA"
            except IndexError:  # This sample isn't in the annotation file
                bcl6trans = "NA"

            # Write the revised sample annotation entry
            outline = [sampleID, copynum, bcl2trans, bcl6trans]
            o.write("\t".join(outline))
            o.write(os.linesep)


# Since we don't have any SV info, just symlink sample annotation file
rule _lymphgen_add_sv_blank:
    input:
        sample_annotation = str(rules._lymphgen_input_cnv.output.sample_annotation)
    output:
        sample_annotation = CFG["dirs"]["lymphgen_input"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sample_annotation.{cnvs_wc}.{sv_wc}.tsv"
    group:
        "lymphgen"
    wildcard_constraints:
        sv_wc = "no_sv"
    run:
        op.relative_symlink(input.sample_annotation, output.sample_annotation, in_module = True)    
        
# STEP 5: RUN LYMPHGEN
 # With CNVs, with A53
rule _lymphgen_run_cnv_A53:
    input:
        sample_annotation = str(rules._lymphgen_add_sv.output.sample_annotation),
        mutation_flat = str(rules._lymphgen_input_cnv.output.mutation_flat),
        gene_list = str(rules._lymphgen_input_cnv.output.gene_list),
        cnv_flat = str(rules._lymphgen_input_cnv.output.cnv_flat),
        cnv_arm = str(rules._lymphgen_input_cnv.output.cnv_arm)
    output:
        result = CFG["dirs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.results.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    log:
        stderr = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stderr.log",
        stdout = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stdout.log"
    params:
        lymphgen_path = CFG["options"]["lymphgen_run"]["lymphgen_path"]
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['optparse']
    wildcard_constraints:
        cnvs_wc = "with_cnvs",
        A53_wc = "with_A53"
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list} -c {input.cnv_flat}
        -a {input.cnv_arm} -o {output.result} > {log.stdout} 2> {log.stderr} """)

# With CNVs, no A53
rule _lymphgen_run_cnv_noA53:
    input:
        sample_annotation = str(rules._lymphgen_add_sv.output.sample_annotation),
        mutation_flat = str(rules._lymphgen_input_cnv.output.mutation_flat),
        gene_list = str(rules._lymphgen_input_cnv.output.gene_list),
        cnv_flat = str(rules._lymphgen_input_cnv.output.cnv_flat),
        cnv_arm = str(rules._lymphgen_input_cnv.output.cnv_arm)
    output:
        result = CFG["dirs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.results.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    log:
        stderr = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stderr.log",
        stdout = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stdout.log"
    params:
        lymphgen_path = CFG["options"]["lymphgen_run"]["lymphgen_path"]
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['optparse']
    wildcard_constraints:
        cnvs_wc = "with_cnvs",
        A53_wc = "no_A53"
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list} -c {input.cnv_flat}
        -a {input.cnv_arm} -o {output.result} --no_A53 > {log.stdout} 2> {log.stderr} """)

# No CNVs
rule _lymphgen_run_no_cnv:
    input:
        sample_annotation = str(rules._lymphgen_add_sv.output.sample_annotation),
        mutation_flat = str(rules._lymphgen_input_no_cnv.output.mutation_flat),
        gene_list = str(rules._lymphgen_input_cnv.output.gene_list)
    output:
        result = CFG["dirs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.results.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    log:
        stderr = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stderr.log",
        stdout = CFG["logs"]["lymphgen_run"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen.{cnvs_wc}.{sv_wc}.{A53_wc}.stdout.log"
    params:
        lymphgen_path = CFG["options"]["lymphgen_run"]["lymphgen_path"]
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['optparse']
    wildcard_constraints:
        cnvs_wc = "no_cnvs"
    shell:
        op.as_one_line("""
        Rscript --vanilla {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list}
        -o {output.result} > {log.stdout} 2> {log.stderr}""")

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lymphgen_output_txt:
    input:
        txt = str(rules._lymphgen_run_cnv_A53.output.result)
    output:
        txt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lymphgen_calls.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    group:
        "lymphgen"
    run:
        op.relative_symlink(input.txt, output.txt, in_module = True)
        
rule _lymphgen_empty_output: 
    output: 
        txt = CFG["dirs"]["outputs"]


# Generates the target sentinels for each run, which generate the symlinks

# Set the applicable wildcards, based on the provided input files
# Are we running LymphGen with CNV/SV data?

# runs with CNV if the cnv_path in the RUNS table is not empty
# runs with SV if there is an sv_path or FISH data for the sample (previous functions alter 
# these values based on whether the sample is genome or capture)

RUNS_CNV_SV = RUNS[((RUNS['has_sv'] == True) | (RUNS['has_fish'] == True)) & (RUNS['cnv_path'] != CFG["inputs"]["sample_seg"]["empty"])]

RUNS_NO_CNV_SV = RUNS[((RUNS['has_sv'] == True) | (RUNS['has_fish'] == True)) & (RUNS['cnv_path'] == CFG["inputs"]["sample_seg"]["empty"])]

RUNS_CNV_NO_SV = RUNS[((RUNS['has_sv'] == False) & (RUNS['has_fish'] == False)) & (RUNS['cnv_path'] != CFG["inputs"]["sample_seg"]["empty"])]

RUNS_NO_CNV_NO_SV = RUNS[((RUNS['has_sv'] == False) & (RUNS['has_fish'] == False)) & (RUNS['cnv_path'] == CFG["inputs"]["sample_seg"]["empty"])]


rule _lymphgen_all:
    input:
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_CNV_SV["tumour_seq_type"],
                genome_build=RUNS_CNV_SV["tumour_genome_build"],
                tumour_id=RUNS_CNV_SV["tumour_sample_id"],
                normal_id=RUNS_CNV_SV["normal_sample_id"],
                pair_status=RUNS_CNV_SV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["with_cnvs", "no_cnvs"],
            sv_wc = ["with_sv", "no_sv"],
            A53_wc = ["with_A53", "no_A53"],
        ),
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_NO_CNV_SV["tumour_seq_type"],
                genome_build=RUNS_NO_CNV_SV["tumour_genome_build"],
                tumour_id=RUNS_NO_CNV_SV["tumour_sample_id"],
                normal_id=RUNS_NO_CNV_SV["normal_sample_id"],
                pair_status=RUNS_NO_CNV_SV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["no_cnvs"],
            sv_wc = ["with_sv", "no_sv"],
            A53_wc = ["no_A53"],
        ),
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_CNV_NO_SV["tumour_seq_type"],
                genome_build=RUNS_CNV_NO_SV["tumour_genome_build"],
                tumour_id=RUNS_CNV_NO_SV["tumour_sample_id"],
                normal_id=RUNS_CNV_NO_SV["normal_sample_id"],
                pair_status=RUNS_CNV_NO_SV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["with_cnvs", "no_cnvs"],
            sv_wc = ["no_sv"],
            A53_wc = ["with_A53", "no_A53"],
        ),
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_NO_CNV_NO_SV["tumour_seq_type"],
                genome_build=RUNS_NO_CNV_NO_SV["tumour_genome_build"],
                tumour_id=RUNS_NO_CNV_NO_SV["tumour_sample_id"],
                normal_id=RUNS_NO_CNV_NO_SV["normal_sample_id"],
                pair_status=RUNS_NO_CNV_NO_SV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["no_cnvs"],
            sv_wc = ["no_sv"],
            A53_wc = ["no_A53"],
        )



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
