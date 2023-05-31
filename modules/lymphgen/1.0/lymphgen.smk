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

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["lymphgen"]`
CFG = op.setup_module(
    name = "lymphgen",
    version = "1.0",
    subdirectories = ["inputs", "reformat_seg", "lymphgen_input", "add_svs", "lymphgen_run", "composite_other", "outputs"],
)

def _find_best_seg(wildcards):
    this_tumor = op.filter_samples(RUNS, genome_build = wildcards.genome_build, tumour_sample_id = wildcards.tumour_id, normal_sample_id = wildcards.normal_id, pair_status = wildcards.pair_status, seq_type = wildcards.seq_type)
    return this_tumor.cnv_path

def _find_best_svar_master(wildcards):
    this_tumor = op.filter_samples(RUNS, genome_build = wildcards.genome_build, tumour_sample_id = wildcards.tumour_id, normal_sample_id = wildcards.normal_id, pair_status = wildcards.pair_status, seq_type = wildcards.seq_type)
    return this_tumor.svar_master_path

def best_seg(mod_config):
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

    all = mod_config["runs"]
    return(all)


RUNS = best_seg(CFG)


# Define rules to be run locally when using a compute cluster
# I put everything under this, since these rules don't take very long
localrules:
    _install_lgenic,
    _lymphgen_input_cnv,
    _lymphgen_input_no_cnv,
    _lymphgen_add_sv,
    _lymphgen_add_sv_blank,
    _lymphgen_run_cnv_A53,
    _lymphgen_run_cnv_noA53, 
    _lymphgen_run_no_cnv,
    _lymphgen_reformat_seg,
    _lymphgen_flag_comp,
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
        download_url=$(curl --silent "https://api.github.com/repos/ckrushton/LGenIC/releases/latest" | grep 'tarball_url' | sed 's/.*:[ ]//' | sed 's/,$//' | sed 's/"//g');
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

# With CNVs
rule _lymphgen_input_cnv:
    input:
        maf = str(rules._lymphgen_input_maf.output.maf),
        seg = str(rules._lymphgen_reformat_seg.output.seg),
        # Software and resource dependencies from LGenIC
        lgenic_script = str(rules._install_lgenic.output.lgenic_script),
        lymphgen_genes = str(rules._install_lgenic.output.lymphgen_genes),
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
        seq_type = CFG["options"]["lymphgen_input"]["seq_type"],
        outprefix = "{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}",
        logratio = "--log2" if CFG["options"]["lymphgen_input"]["use_log_ratio"].lower() == "true" else ""
    conda:
        CFG['conda_envs']['sorted_containers']
    wildcard_constraints:
        cnvs_wc = "with_cnvs"
    shell:
        op.as_one_line("""
        python {input.lgenic_script} --lymphgen_genes {input.lymphgen_genes} --sequencing_type {params.seq_type} --outdir $(dirname {output.sample_annotation})
        --outprefix {params.outprefix} -v INFO --maf {input.maf} --entrez_ids {input.hugo2entrez} --cnvs {input.seg} {params.logratio} --genes {input.gene_coords} --arms {input.arm_coords}
        > {log.stdout} 2> {log.stderr}
        """)

# No CNVs
rule _lymphgen_input_no_cnv:
    input:
        maf = str(rules._lymphgen_input_maf.output.maf),
        # Software and resource dependencies from LGenIC
        lgenic_script = str(rules._install_lgenic.output.lgenic_script),
        lymphgen_genes = str(rules._install_lgenic.output.lymphgen_genes),
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
        seq_type = CFG["options"]["lymphgen_input"]["seq_type"],
        outprefix = "{tumour_id}--{normal_id}--{pair_status}.{cnvs_wc}"
    conda:
        CFG['conda_envs']['sorted_containers']
    wildcard_constraints:
        cnvs_wc = "no_cnvs"
    shell:
        op.as_one_line("""
        python {input.lgenic_script} --lymphgen_genes {input.lymphgen_genes} --sequencing_type {params.seq_type} --outdir $(dirname {output.sample_annotation})
        --outprefix {params.outprefix} -v INFO --maf {input.maf} --entrez_ids {input.hugo2entrez} > {log.stdout} 2> {log.stderr}
        """)

# STEP 4: Add SV information (if availible)

# Obtain the path to the GAMBLR conda environment
md5hash = hashlib.md5()
if workflow.conda_prefix:
    conda_prefix = workflow.conda_prefix
else:
    conda_prefix = os.path.abspath(".snakemake/conda")

md5hash.update(conda_prefix.encode())
f = open(config["lcr-modules"]["lymphgen"]["conda_envs"]["gamblr"], 'rb')
md5hash.update(f.read())
f.close()
h = md5hash.hexdigest()
GAMBLR = glob.glob(conda_prefix + "/" + h[:8] + "*")
for file in GAMBLR: 
    if os.path.isdir(file): 
        GAMBLR = file

rule _lymphgen_install_GAMBLR:
    params:
        branch = ", ref = \"" + CFG['inputs']['gamblr_branch'] + "\"" if CFG['inputs']['gamblr_branch'] != "" else "",
        config_url = CFG["inputs"]["gamblr_config_url"]
    output:
        installed = directory(GAMBLR + "/lib/R/library/GAMBLR"),
        config = "gamblr.yaml"
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['gamblr']
    shell:
        op.as_one_line("""
        wget -qO {output.config} {params.config_url} &&
        R -q -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR"{params.branch})'
        """)

#print(CFG["inputs"]["sample_sv_info"]["svar_master"])

def get_svar_master(wildcards):
    possible_output = (expand(config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["svar_master"], zip,
                                        tumour_id=wildcards.tumour_id,
                                        normal_id=wildcards.normal_id,
                                        pair_status=wildcards.pair_status,
                                        genome_build=wildcards.genome_build,
                                        seq_type=wildcards.seq_type))

    
    if os.path.exists(possible_output[0]):
        return possible_output[0]
    else:
        return config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["empty_svar"]


def get_svar(wildcards):
    if config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["svar_master"] == config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["empty_svar"]:
        return config["lcr-modules"]["lymphgen"]["inputs"]["sample_sv_info"]["empty_svar"]
    else:
        return get_svar_master(wildcards)

rule _lymphgen_input_sv:
    input:
        fish = CFG["inputs"]["sample_sv_info"]["fish"],
        svar_master = get_svar, 
        gamblr = ancient(rules._lymphgen_install_GAMBLR.output.installed)
    output:
        sv = CFG["dirs"]["inputs"] + "sv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.tsv"
    group:
        "lymphgen"
    conda:
        CFG['conda_envs']['gamblr']
    script:
        config["lcr-modules"]["lymphgen"]['options']['lymphgen_run']['lymphgen_sv']

rule _lymphgen_add_sv:
    input:
        sample_annotation = str(rules._lymphgen_input_cnv.output.sample_annotation),
        bcl2_bcl6_sv = str(rules._lymphgen_input_sv.output.sv)
    output:
        sample_annotation = CFG["dirs"]["add_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sample_annotation.{cnvs_wc}.{sv_wc}.tsv"
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

        # Add SV info to sample annotation file
        with open(input.sample_annotation) as f, open(output.sample_annotation, "w") as o:
            for line in f:
                line = line.rstrip("\n").rstrip("\r")  # Handle line endings
                if line.startswith("Sample.ID\tCopy.Number"):  # Header line
                    o.write(line)
                    o.write(os.linesep)
                    continue

                cols = line.split("\t")
                try:
                    sampleID = cols[0]
                    copynum = cols[1]
                except IndexError as e:
                    raise AttributeError("Input sample annotaion file \'%s\' appears to be malformed" % input.sample_annotation) from e
                # Find the matching SampleID in the SV annotation file
                sampleEntry = loaded_sv.loc[loaded_sv[params.sampleIDcolname] == sampleID]
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
        sample_annotation = CFG["dirs"]["add_svs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}_sample_annotation.{cnvs_wc}.{sv_wc}.tsv"
    group:
        "lymphgen"
    wildcard_constraints:
        sv_wc = "no_sv"
    run:
        op.relative_symlink(input.sample_annotation, output.sample_annotation)


# STEP 5: RUN LYMPHGEN

def _get_sample_annotation(wildcards):
    if wildcards.sv_wc == "has_sv":
        return str(rules._lymphgen_add_sv.output.sample_annotation)
    else:
        return str(rules._lymphgen_add_sv_blank.output.sample_annotation)

# With CNVs, with A53
rule _lymphgen_run_cnv_A53:
    input:
        sample_annotation = _get_sample_annotation,
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
        Rscript {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list} -c {input.cnv_flat}
        -a {input.cnv_arm} -o {output.result} > {log.stdout} 2> {log.stderr} """)

# With CNVs, no A53
rule _lymphgen_run_cnv_noA53:
    input:
        sample_annotation = _get_sample_annotation,
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
        Rscript {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list} -c {input.cnv_flat}
        -a {input.cnv_arm} -o {output.result} --no_A53 > {log.stdout} 2> {log.stderr} """)

# No CNVs
rule _lymphgen_run_no_cnv:
    input:
        sample_annotation = _get_sample_annotation,
        mutation_flat = str(rules._lymphgen_input_no_cnv.output.mutation_flat),
        gene_list = str(rules._lymphgen_input_no_cnv.output.gene_list)
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
        Rscript {params.lymphgen_path} -m {input.mutation_flat} -s {input.sample_annotation} -g {input.gene_list}
        -o {output.result} > {log.stdout} 2> {log.stderr}""")

# STEP 6. Flag samples that fall in the composite "Dead Zone"
rule _lymphgen_flag_comp:
    input:
        txt = str(rules._lymphgen_run_cnv_A53.output.result)
    output:
        txt = CFG["dirs"]["composite_other"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.lymphgen_composite_other.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    group:
        "lymphgen"
    run:
        loaded_calls = pandas.read_csv(input.txt, sep="\t")

        # Subset down to cases not classified by LymphGen (i.e. "Other")
        othercases = loaded_calls.loc[loaded_calls["Subtype.Prediction"] == "Other"]

        # Grab the confidence columns. These will change depending on which subgroups were included in the classification
        confcols = list(x for x in othercases.columns if x.startswith("Confidence"))
        confvalues = othercases[confcols] > 0.5
        compother = othercases[confvalues.any(1)]
        compother.to_csv(output.txt, sep = "\t")


# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _lymphgen_output_txt:
    input:
        txt = str(rules._lymphgen_run_cnv_A53.output.result),
        comp = str(rules._lymphgen_flag_comp.output.txt)
    output:
        txt = CFG["dirs"]["outputs"] + "{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}/lymphgen_calls.{cnvs_wc}.{sv_wc}.{A53_wc}.tsv"
    group:
        "lymphgen"
    run:
        op.relative_symlink(input.txt, output.txt)


# Generates the target sentinels for each run, which generate the symlinks

# Set the applicable wildcards, based on the provided input files
# Are we running LymphGen with CNV/SV data?

if "sample_sv_info" in CFG["inputs"] and CFG["inputs"]["sample_sv_info"] != "" and CFG["inputs"]["sample_sv_info"] != "None":
    sv_wc = ["with_sv", "no_sv"]
else:
    sv_wc = ["no_sv"]

# Are we including A53 subgroup in the output?
if "sample_seg" in CFG["inputs"] and CFG["inputs"]["sample_seg"] != "" and CFG["inputs"]["sample_seg"] != "None":  # We have CNV calls
    a53_wc = ["with_A53", "no_A53"]
else:
    a53_wc = ["no_A53"]

RUNS_CNV = RUNS[RUNS['cnv_path'] != CFG["inputs"]["sample_seg"]["empty"]]

RUNS_NO_CNV = RUNS[RUNS['cnv_path'] == CFG["inputs"]["sample_seg"]["empty"]]

rule _lymphgen_all:
    input:
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_CNV["tumour_seq_type"],
                genome_build=RUNS_CNV["tumour_genome_build"],
                tumour_id=RUNS_CNV["tumour_sample_id"],
                normal_id=RUNS_CNV["normal_sample_id"],
                pair_status=RUNS_CNV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["with_cnvs", "no_cnvs"],
            sv_wc = sv_wc,
            A53_wc = a53_wc,
        ),
        expand(
            expand(
                [
                    str(rules._lymphgen_output_txt.output.txt),
                ],
                zip,
                seq_type=RUNS_NO_CNV["tumour_seq_type"],
                genome_build=RUNS_NO_CNV["tumour_genome_build"],
                tumour_id=RUNS_NO_CNV["tumour_sample_id"],
                normal_id=RUNS_NO_CNV["normal_sample_id"],
                pair_status=RUNS_NO_CNV["pair_status"], 
                allow_missing = True
            ),
            cnvs_wc = ["no_cnvs"],
            sv_wc = sv_wc,
            A53_wc = a53_wc,
        )



##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
