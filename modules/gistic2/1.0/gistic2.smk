#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Sierra Gillis
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
# `CFG` is a shortcut to `config["lcr-modules"]["gistic2"]`
CFG = op.setup_module(
    name = "gistic2",
    version = "1.0",
    subdirectories = ["inputs", "prepare_seg", "markers", "gistic2", "outputs"],
)

# Define rules to be run locally when using a compute cluster
localrules:
    _gistic2_input_seg,
    _gistic2_input_subsetting_categories,
    _gistic2_prepare_seg,
    _gistic2_make_markers,
    _gistic2_download_ref,
    _gistic2_output,
    _gistic2_all


##### RULES #####
launch_date = datetime.today().strftime('%Y-%m')

# Symlinks the input files into the module results directory (under '00-inputs/')
rule _gistic2_input_seg:
    input:
        seg = CFG["inputs"]["seg"]
    output:
        seg = CFG["dirs"]["inputs"] + "{seq_type}--projection/all--{projection}.seg"
    group: 
        "input_and_format"
    run:
        op.absolute_symlink(input.seg, output.seg)

rule _gistic2_input_subsetting_categories:
    input:
        subsetting_categories = CFG["inputs"]["subsetting_categories"]
    output:
        subsetting_categories = CFG["dirs"]["inputs"] + "sample_sets/subsetting_categories.tsv"
    group: 
        "input_and_format"
    run:
        op.absolute_symlink(input.subsetting_categories, output.subsetting_categories)

# Download refgene reference MAT file
rule _gistic2_download_ref:
    output:
        refgene_mat = CFG["dirs"]["inputs"] + "references/{projection}.refgene.mat"
    group: 
        "input_and_format"
    params:
        url = "http://bcgsc.ca/downloads/morinlab/gistic2_references/{projection}.refgene.mat",
        folder = CFG["dirs"]["inputs"] + "references"
    shell:
        op.as_one_line("""
        wget -P {params.folder} {params.url} 
        """)

# Merges capture and genome seg files if available, and subset to the case_set provided
checkpoint _gistic2_prepare_seg:
    input:
        seg = expand(str(rules._gistic2_input_seg.output.seg),
                    allow_missing=True,
                    seq_type=CFG["samples"]["seq_type"].unique()
                    ),
        subsetting_categories = str(rules._gistic2_input_subsetting_categories.output.subsetting_categories)
    output:
        directory(CFG["dirs"]["prepare_seg"] + "{case_set}--{projection}--{launch_date}")
    log:
        log = CFG["logs"]["prepare_seg"] + "{case_set}--{projection}--{launch_date}.log"
    group: 
        "input_and_format"
    params:
        launch_date = launch_date,
        case_set = CFG["case_set"],
        seq_type = CFG["samples"]["seq_type"].unique(),
        metadata = CFG["samples"][["sample_id","seq_type","genome_build","cohort","pathology","unix_group","time_point"]].to_numpy(na_value='')
    script:
        config["lcr-modules"]["gistic2"]["prepare_seg"]

def _gistic2_get_seg_with_md5sum(wildcards):
    CFG = config["lcr-modules"]["gistic2"]
    checkpoint_output = checkpoints._gistic2_prepare_seg.get(**wildcards).output[0]
    print("checkpoint_output", checkpoint_output)
    sums, = glob_wildcards(checkpoint_output + "/{md5sum}.seg")
    print("sums", sums)
    segs = expand(
            CFG["dirs"]["prepare_seg"] + "{case_set}--{projection}--{launch_date}/{md5sum}.seg", md5sum = sums, launch_date=launch_date, allow_missing=True)
    
    return(segs)

# Create a markers file that has every segment start and end that appears in the seg file
rule _gistic2_make_markers:
    input:
        seg = _gistic2_get_seg_with_md5sum
    output:
        temp_markers = temp(CFG["dirs"]["markers"] + "temp_markers--{case_set}--{projection}--{launch_date}--{md5sum}.txt"),
        markers = CFG["dirs"]["markers"] + "markers--{case_set}--{projection}--{launch_date}--{md5sum}.txt"
    log:
        stderr = CFG["logs"]["markers"] + "{case_set}--{projection}--{launch_date}--{md5sum}--markers.stderr.log"
    shell: 
        op.as_one_line("""
        sed '1d' {input.seg} | cut -f2,3 > {output.temp_markers}
            &&
        sed '1d' {input.seg} | cut -f2,4 >> {output.temp_markers}
            &&
        sort -V -k1,1 -k2,2nr {output.temp_markers} | uniq | nl > {output.markers} 
        2> {log.stderr}
        """)

# Run gistic2 for a single seq_type (capture, genome) for the confidence thresholds listed in the config
rule _gistic2_run:
    input:
        seg = _gistic2_get_seg_with_md5sum,
        refgene_mat = str(rules._gistic2_download_ref.output.refgene_mat),
        markers = str(rules._gistic2_make_markers.output.markers)
    output:
        all_data_by_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/all_lesions.conf_{conf}.txt",
        amp_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/amp_genes.conf_{conf}.txt",
        del_genes = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/del_genes.conf_{conf}.txt",
        scores = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/scores.gistic"
    log:
        stdout = CFG["logs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/gistic2.stdout.log",
        stderr = CFG["logs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/gistic2.stderr.log"
    params:
        base_dir = CFG["dirs"]["gistic2"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/",
        opts = CFG["options"]["gistic2_run"]
    conda:
        CFG["conda_envs"]["gistic2"]
    threads:
        CFG["threads"]["gistic2_run"]
    resources:
        **CFG["resources"]["gistic2_run"]
    shell:
        op.as_one_line("""
        gistic2 -b {params.base_dir} -seg {input.seg} -mk {input.markers} -refgene {input.refgene_mat} 
        -genegistic 1 -broad 1 -savegene 1 -conf 0.{wildcards.conf} -v 30 -saveseg 0 -savedata 0
        {params.opts}
        > {log.stdout} 2> {log.stderr}
        """)

# Symlinks the final output files into the module results directory (under '99-outputs/')
rule _gistic2_output:
    input:
        all_data_by_genes = str(rules._gistic2_run.output.all_data_by_genes),
        all_lesions = str(rules._gistic2_run.output.all_lesions),
        amp_genes = str(rules._gistic2_run.output.amp_genes),
        del_genes = str(rules._gistic2_run.output.del_genes),
        scores = str(rules._gistic2_run.output.scores)
    output:
        all_data_by_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/all_data_by_genes.txt",
        all_lesions = CFG["dirs"]["outputs"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/all_lesions.conf_{conf}.txt",
        amp_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/amp_genes.conf_{conf}.txt",
        del_genes = CFG["dirs"]["outputs"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/del_genes.conf_{conf}.txt",
        scores = CFG["dirs"]["outputs"] + "{case_set}--{projection}/{launch_date}--{md5sum}/conf_{conf}/scores.gistic"

    run:
        op.relative_symlink(input.all_data_by_genes, output.all_data_by_genes, in_module= True),
        op.relative_symlink(input.all_lesions, output.all_lesions, in_module= True),
        op.relative_symlink(input.amp_genes, output.amp_genes, in_module= True),
        op.relative_symlink(input.del_genes, output.del_genes, in_module= True),
        op.relative_symlink(input.scores, output.scores, in_module= True)

def _gistic2_get_md5sums(wildcards):
    checkpoint_output = checkpoints._gistic2_prepare_seg.get(**wildcards).output[0]
    sums, = glob_wildcards(checkpoint_output + launch_date + "--{md5sum}.seg")

    return(sums)

rule _gistic2_all:
    input:
        expand(
            [
                str(rules._gistic2_output.output.all_data_by_genes),
                str(rules._gistic2_output.output.all_lesions),
                str(rules._gistic2_output.output.amp_genes),
                str(rules._gistic2_output.output.del_genes),
                str(rules._gistic2_output.output.scores)
            ],              
            conf = CFG["options"]["conf_level"],
            projection = CFG["projections"],
            case_set = CFG["case_set"],
            launch_date = launch_date,
            md5sum = _gistic2_get_md5sums
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
