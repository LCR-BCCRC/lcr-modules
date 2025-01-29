#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Laura Hilton
# Module Author:    Laura Hilton
# Contributors:     Sierra Gillis


##### SETUP #####

# Import package with useful functions for developing analysis modules
import sys, os
from os.path import join
import glob
import oncopipe as op
from datetime import datetime

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
# `CFG` is a shortcut to `config["lcr-modules"]["dlbclass"]`
CFG = op.setup_module(
    name = "dlbclass",
    version = "1.0",
    subdirectories = ["inputs", "prepare_inputs", "dlbclass", "outputs"],
)

# All rules should be run locally
localrules:
    _dlbclass_input_maf,
    _dlbclass_input_seg,
    _dlbclass_input_sv,
    _dlbclass_input_subsetting_categories,
    _dlbclass_download_dlbclass,
    _dlbclass_download_refs,
    _dlbclass_maf_to_gsm, 
    _dlbclass_seg_to_gsm, 
    _dlbclass_sv_to_gsm,
    _dlbclass_combine_gsm,
    _dlbclass_run, 
    _dlbclass_output,
    _dlbclass_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

# Interpret the absolute path to this script so it doesn't get interpreted relative to the module snakefile later.
PREPARE_MAFS =  os.path.abspath(config["lcr-modules"]["dlbclass"]["scripts"]["prepare_mafs"])
PREPARE_SVS =  os.path.abspath(config["lcr-modules"]["dlbclass"]["scripts"]["prepare_svs"])

# Get date as used in DLBCLass output files
TODAY = datetime.now().strftime("%d%b%Y")

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _dlbclass_input_maf:
    input:
        maf = CFG["inputs"]["master_maf"]
    output:
        maf = CFG["dirs"]["inputs"] + "variants/{sample_set}--{launch_date}/{seq_type}.input.maf"
    run:
        op.absolute_symlink(input.maf, output.maf)
        
rule _dlbclass_input_seg:
    input:
        seg = CFG["inputs"]["master_seg"] if CFG["inputs"]["master_seg"] != "" else CFG["inputs"]["empty_seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "variants/{sample_set}--{launch_date}/{seq_type}.input.seg", 
    run:
        op.absolute_symlink(input.seg, output.seg)
        
rule _dlbclass_input_sv:
    input:
        sv = CFG["inputs"]["master_sv"] if CFG["inputs"]["master_sv"] != "" else CFG["inputs"]["empty_sv"]
    output:
        sv = CFG["dirs"]["inputs"] + "variants/{sample_set}--{launch_date}/{seq_type}.input.sv"
    run:
        op.absolute_symlink(input.sv, output.sv)

# Symlinks the subsetting categories input file into the module results directory (under '00-inputs/')
rule _dlbclass_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

# Download the dlbclass
rule _dlbclass_download_dlbclass:
    output:
        dlbclass_compressed = temp(CFG["dirs"]["inputs"] + "DLBclass-tool.tar.gz"),
        dlbclass = CFG["dirs"]["inputs"] + "DLBclass-tool/dlbclass.py", 
        seg2gsm = CFG["dirs"]["inputs"] + "DLBclass-tool/src/seg2gsm.py", 
        maf2gsm = CFG["dirs"]["inputs"] + "DLBclass-tool/src/maf2gsm.py", 
        combine2gsm = CFG["dirs"]["inputs"] + "DLBclass-tool/src/combine2gsm.py",
        feature_order = CFG["dirs"]["inputs"] + "DLBclass-tool/gsm/feature_order.19Aug2024.txt",
        pub_gsm = CFG["dirs"]["inputs"] + "DLBclass-tool/gsm/DLBclass_published.GSM.tsv"
    params: 
        dlbclass_release = CFG["inputs"]["dlbclass_release"]
    conda:
        CFG["conda_envs"]["wget"]
    shell:
        op.as_one_line("""
        wget -cO - {params.dlbclass_release} > {output.dlbclass_compressed} && tar -C $(dirname {output.dlbclass}) -xf {output.dlbclass_compressed};
        mv $(dirname {output.dlbclass})/DLBclass-tool-*/* $(dirname {output.dlbclass})/ && rm -r $(dirname {output.dlbclass})/DLBclass-tool-*;
        """)

# Download the reference files
rule _dlbclass_download_refs:
    params: 
        focal_cnv = CFG["inputs"]["focal_cnv"],
        arm_cnv = CFG["inputs"]["arm_cnv"],
        blacklist_cnv = CFG["inputs"]["blacklist_cnv"]
    output:
        focal_cnv = CFG["dirs"]["inputs"] + "refs/" + CFG["inputs"]["focal_cnv"].split("/")[-1],
        arm_cnv = CFG["dirs"]["inputs"] + "refs/" + CFG["inputs"]["arm_cnv"].split("/")[-1],
        blacklist_cnv = CFG["dirs"]["inputs"] + "refs/" + CFG["inputs"]["blacklist_cnv"].split("/")[-1]
    conda:
        CFG["conda_envs"]["wget"] 
    shell:
        op.as_one_line("""
        wget
        -O - {params.focal_cnv} | 
        grep -v "19q13.32_2:DEL" | 
        sed 's/19q13.32_1:DEL/19q13.32:DEL/g' > {output.focal_cnv}
            &&
        wget
        -O {output.arm_cnv}
        {params.arm_cnv}
            &&
        wget
        -O {output.blacklist_cnv}
        {params.blacklist_cnv}
        """)


# Prepare the maf file for the input to dlbclass
checkpoint _dlbclass_prepare_maf:
    input:
        maf = expand(
                    str(rules._dlbclass_input_maf.output.maf),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._dlbclass_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/maf.done"
    log:
        CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/prepare_maf.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "dlbclass",
        metadata_cols = CFG["samples"],
        metadata_dim = CFG["samples"].shape,
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

# Prepare the seg file for the input to dlbclass        
checkpoint _dlbclass_prepare_seg:
    input:
        seg = expand(
                    str(rules._dlbclass_input_seg.output.seg),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._dlbclass_input_subsetting_categories.output.subsetting_categories)
    output:
        CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/seg.done"
    log:
        CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/prepare_seg.log"
    conda:
        CFG["conda_envs"]["prepare_mafs"]
    params:
        include_non_coding = str(CFG["include_non_coding"]).upper(),
        mode = "dlbclass",
        metadata_cols = CFG["samples"],
        metadata_dim = CFG["samples"].shape,
        metadata = CFG["samples"].to_numpy(na_value='')
    script:
        PREPARE_MAFS

# Prepare MAF GSM file
def _get_input_maf(wildcards):
    CFG = config["lcr-modules"]["dlbclass"]
    checkpoint_output = os.path.dirname(str(checkpoints._dlbclass_prepare_maf.get(**wildcards).output[0]))
    SUMS, = glob_wildcards(checkpoint_output +"/{md5sum_maf}.maf")
    inputs = expand(
        [
            CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{md5sum_maf}.maf",
            CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{md5sum_maf}.maf.content"
        ],
        md5sum_maf = SUMS, 
        allow_missing = True
        )
    return inputs

rule _dlbclass_maf_to_gsm:
    input:
        _get_input_maf, 
        script = str(rules._dlbclass_download_dlbclass.output.maf2gsm)
    output:
        maf_gsm = CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{sample_set}." + TODAY + ".MAF.GSM.tsv"
    log: CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/maf2gsm.log"
    conda:
        CFG["conda_envs"]["dlbclass"]
    shell:
        op.as_one_line("""
        python3 {input.script}   
            --id {wildcards.sample_set} 
            -s {input[1]} 
            -m {input[0]}
            -o $(dirname {output.maf_gsm}) 
            > {log} 2>&1 
        """)

# Prepare CNV GSM file
def _get_input_seg(wildcards):
    CFG = config["lcr-modules"]["dlbclass"]
    seg_output = os.path.dirname(str(checkpoints._dlbclass_prepare_seg.get(**wildcards).output[0]))
    SEG_SUMS, = glob_wildcards(seg_output +"/{md5sum_seg}.seg")
    maf_output = os.path.dirname(str(checkpoints._dlbclass_prepare_maf.get(**wildcards).output[0]))
    MAF_SUMS, = glob_wildcards(maf_output +"/{md5sum_maf}.maf")
    inputs = expand(
        [
            CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{md5sum_seg}.seg",
            CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{md5sum_maf}.maf.content"
        ],
        md5sum_seg = SEG_SUMS, 
        md5sum_maf = MAF_SUMS,
        allow_missing = True
        )
    return inputs

rule _dlbclass_seg_to_gsm:
    input:
        _get_input_seg,
        _get_input_maf,
        arm_cnv = str(rules._dlbclass_download_refs.output.arm_cnv), 
        focal_cnv = str(rules._dlbclass_download_refs.output.focal_cnv), 
        blacklist_cnv = str(rules._dlbclass_download_refs.output.blacklist_cnv),
        script = str(rules._dlbclass_download_dlbclass.output.seg2gsm)
    output:
        cnv_gsm = CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{sample_set}." + TODAY + ".CNV.GSM.tsv"
    log: CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/seg2gsm.log"
    conda:
        CFG["conda_envs"]["dlbclass"]
    shell:
        op.as_one_line("""
        python3 {input.script} 
            -i {wildcards.sample_set} 
            -s {input[3]} 
            -v {input[0]}
            -x {input.blacklist_cnv}
            -a {input.arm_cnv}
            -f {input.focal_cnv}
            -o $(dirname {output.cnv_gsm}) 
            -g hg19
            > {log} 2>&1
        """)

rule _dlbclass_sv_to_gsm:
    input:
        _get_input_maf, 
        sv = expand(
                str(rules._dlbclass_input_sv.output.sv),
                allow_missing=True,
                seq_type=CFG["samples"]["seq_type"].unique()
                )       
    output:
        sv_gsm = CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{sample_set}." + TODAY + ".SV.GSM.tsv"
    params: 
        pos_value = CFG["options"]["sv2gsm"]["POS_value"],
        neg_value = CFG["options"]["sv2gsm"]["NEG_value"]
    log: CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/sv2gsm.log"
    conda:
        CFG["conda_envs"]["dlbclass"]
    script:
        PREPARE_SVS
        

def get_combine_param(wildcards):
    CFG = config["lcr-modules"]["dlbclass"]
    if(wildcards.sv_wc == "with_sv"):
        sv_param = "-v " + str(expand(rules._dlbclass_sv_to_gsm.output.sv_gsm, **wildcards, allow_missing = True)[0])
    else:
        sv_param = ""
        
    if(wildcards.cnv_wc == "with_cnv"):
        cnv_param = "-c " + str(expand(rules._dlbclass_seg_to_gsm.output.cnv_gsm, **wildcards, allow_missing = True)[0])
    else:
        cnv_param = ""
        
    return cnv_param + " " + sv_param


# def get_combine_input(wildcards):
#     CFG = config["lcr-modules"]["dlbclass"]
#     if(wildcards.sv_wc == "with_sv"):
#         sv_param = "-v " + str(rules._dlbclass_sv_to_gsm.output.sv_gsm)
#     else:
#         sv_param = ""
        
#     if(wildcards.cnv_wc == "with_cnv"):
#         cnv_param = "-c " + str(rules._dlbclass_seg_to_gsm.output.cnv_gsm)
#     else:
#         cnv_param = ""
        
#     return cnv_param + " " + sv_param

# Actual dlbclass run
rule _dlbclass_combine_gsm:
    input:
        sv = lambda w: str(rules._dlbclass_sv_to_gsm.output.sv_gsm) if w.sv_wc == "with_sv" else config["lcr-modules"]["dlbclass"]["inputs"]["empty_sv"],
        seg = lambda w: str(rules._dlbclass_seg_to_gsm.output.cnv_gsm) if w.cnv_wc == "with_cnv" else config["lcr-modules"]["dlbclass"]["inputs"]["empty_seg"],
        maf_gsm = str(rules._dlbclass_maf_to_gsm.output.maf_gsm), 
        feature_order = str(rules._dlbclass_download_dlbclass.output.feature_order), 
        script = str(rules._dlbclass_download_dlbclass.output.combine2gsm), 
        pub_gsm = str(rules._dlbclass_download_dlbclass.output.pub_gsm)
    output:
        gsm = CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/{sample_set}--{sv_wc}--{cnv_wc}." + TODAY + ".GSM.tsv"
    params: 
        get_combine_param
    log: CFG["logs"]["prepare_inputs"] + "{sample_set}--{launch_date}/combine2gsm.{sv_wc}--{cnv_wc}.log"
    conda:
        CFG["conda_envs"]["dlbclass"]
    threads:
        CFG["threads"]["dlbclass"]
    resources:
        **CFG["resources"]["dlbclass"]
    shell:
        op.as_one_line("""
        python3 {input.script} 
        -i {wildcards.sample_set}--{wildcards.sv_wc}--{wildcards.cnv_wc} 
        -m {input.maf_gsm} 
        -f {input.feature_order}
        -p {input.pub_gsm}
        {params[0]}
        -o $(dirname {output.gsm})
        > {log} 2>&1
        """)

rule _dlbclass_run:
    input:
        gsm = str(rules._dlbclass_combine_gsm.output.gsm), 
        script = str(rules._dlbclass_download_dlbclass.output.dlbclass)
    output:
        dlbclass = CFG["dirs"]["dlbclass"] + "{sample_set}--{launch_date}/{sample_set}--{sv_wc}--{cnv_wc}_classified_samples.tsv",
    conda:
        CFG["conda_envs"]["dlbclass"]
    threads:
        CFG["threads"]["dlbclass"]
    resources:
        **CFG["resources"]["dlbclass"]
    shell:
        op.as_one_line("""
        python3 {input.script}
        -i {wildcards.sample_set}--{wildcards.sv_wc}--{wildcards.cnv_wc} 
        -g {input.gsm} 
        -o $(dirname {output.dlbclass})
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _dlbclass_output:
    input:
        tsv = str(rules._dlbclass_run.output.dlbclass)
    output:
        tsv = CFG["dirs"]["outputs"] + "txt/{sample_set}--{launch_date}/{sample_set}--{sv_wc}--{cnv_wc}.dlbclass.tsv"
    run:
        op.relative_symlink(input.tsv, output.tsv, in_module= True)

cnv_wc = ["with_cnv", "no_cnv"] if CFG["inputs"]["master_seg"] != "" else ["no_cnv"]
sv_wc = ["with_sv", "no_sv"] if CFG["inputs"]["master_sv"] != "" else ["no_sv"]

# Generates the target sentinels for each run, which generate the symlinks
rule _dlbclass_all:
    input:
        expand(
            [
                CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/maf.done",
                str(rules._dlbclass_output.output.tsv),
            ],
            sample_set=CFG["sample_set"],
            cnv_wc = cnv_wc, 
            sv_wc = sv_wc,
            launch_date = launch_date)


##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
