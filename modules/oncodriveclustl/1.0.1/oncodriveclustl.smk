#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  N/A
# Module Author:    Manuela Cruz
# Contributors:     N/A


import oncopipe as op
import os
import pandas as pd
import hashlib

### stuff that is not supposed to be here:
configfile: "src/lcr-modules/modules/oncodriveclustl/1.0.1/config/default.yaml"

config["tool_names"] = "oncodriveclustl"
config["pipeline_name"] = "oncodriveclustl"

include: "header_2.1.smk"

config["lcr-modules"]["oncodriveclustl"]["dirs"]["_parent"] = config['lcr-modules']['_shared']['root_output_dir'] + "oncodriveclustl_module_test/"

#### FOR CHECKING MODULE RUNS:
RUN_SAMPLES = PIPELINE_SAMPLES["SAMPLES"]["oncodriveclustl"]
#### add this to code below? kinda necessary to run?
####

if isinstance(config["lcr-modules"]["oncodriveclustl"]["options"]["pathology"],list):
    if isinstance(config["lcr-modules"]["oncodriveclustl"]["options"]["pathology"][0], list):
        run_pathology = sum(config["lcr-modules"]["oncodriveclustl"]["options"]["pathology"], [])
        run_pathology = list(set(run_pathology))
    else:
        run_pathology = config["lcr-modules"]["oncodriveclustl"]["options"]["pathology"]

if config["lcr-modules"]["oncodriveclustl"]["options"]["pathology"] is not None:
    RUN_SAMPLES = op.filter_samples(RUN_SAMPLES, seq_type=config["lcr-modules"]["oncodriveclustl"]["options"]["seq_type"], pathology=run_pathology)
    RUN_IDS = RUN_SAMPLES.patient_id.tolist()
    PATHOLOGY = set(RUN_SAMPLES.pathology.tolist())
    SEQ_TYPES = set(RUN_SAMPLES.seq_type.tolist())

    print(RUN_IDS)
    print(PATHOLOGY)
    print(SEQ_TYPES)

elif config["lcr-modules"]["oncodriveclustl"]["options"]["sample_ids"] is not None:
    list_sample_ids = config["lcr-modules"]["oncodriveclustl"]["options"]["sample_ids"].split(",")
    RUN_SAMPLES = op.filter_samples(RUN_SAMPLES, seq_type=config["lcr-modules"]["oncodriveclustl"]["options"]["seq_type"], sample_id=list_sample_ids)
    RUN_IDS = RUN_SAMPLES.patient_id.tolist()
    PATHOLOGY = set(RUN_SAMPLES.pathology.tolist())
    SEQ_TYPES = set(RUN_SAMPLES.seq_type.tolist())

    print(RUN_IDS)
    print(PATHOLOGY)
    print(SEQ_TYPES)

config["lcr-modules"]["oncodriveclustl"]["samples"] = RUN_SAMPLES


###################################

##### SETUP #####


# Import package with useful functions for developing analysis modules
#import oncopipe as op

# Setup module and store module-specific configuration in `CFG`
# `CFG` is a shortcut to `config["lcr-modules"]["mixcr"]`
CFG = op.setup_module(
    name = "oncodriveclustl",
    version = "1.0",
    subdirectories = ["inputs", "oncodriveclustl", "outputs"],
)

assert type(CFG['coding_sequence'])==bool, (
    "Config value for 'coding_sequence' must be set to a Boolean value (True/False)"
    "Default is True, only CDS regions will be considered as input for CLUSTL."
)

assert CFG["options"]["pathology"] is not None or CFG["options"]["sample_ids"] is not None or CFG["options"]["sample_subset"] is not None, (
    "Please specify pathologies, sample_ids or subgroups to run."
)

assert CFG["options"]["mode"] is not None, (
    "Please specify mode for preparing oncodrive input, choose from: 'pathology', 'samples', 'subset', 'maf'"
)
##### One way to check to see if CONTENTS file exists or can be updated for pathologies


BUILD = CFG["options"]["genome_build"]
contents_dir = CFG["dirs"]["oncodriveclustl"] + f"{BUILD}/contents/"
RUN_SAMPLES = config["lcr-modules"]["oncodriveclustl"]["samples"]


if CFG["options"]["pathology"] is not None:
    CFG["options"]["mode"] = "pathology"

    SEQUENCIES = CFG["options"]["seq_type"]
    PATHOLOGIES = CFG["options"]["pathology"]
    ## Just for testing, ideally pathology input would be in list format
    if not isinstance(PATHOLOGIES,list):
        PATHOLOGIES = PATHOLOGIES.split(' ')
    if not isinstance(PATHOLOGIES[0],list):
        PATHOLOGIES = [PATHOLOGIES]

    FINAL_PATHOLOGIES = {}
    DESIRED_PATHOLOGIES = {}

    print(PATHOLOGIES)

    for x in PATHOLOGIES:
        alphabetical_l = "--".join(sorted(x))
        l = alphabetical_l
        print(l)
        new = True
        if os.path.exists(contents_dir):
            for f in os.listdir(contents_dir):
                if alphabetical_l in f and f.endswith(".contents"):
                    new = False
        if new:
            FINAL_PATHOLOGIES[l] = 1
            DESIRED_PATHOLOGIES[l] = 1 
        if not new:
            version_count = 0
            if os.path.exists(contents_dir):
                for f in os.listdir(contents_dir):
                    file_prefix = f"{alphabetical_l}_v"
                    if file_prefix in f and f.endswith(".contents"):
                        version_count += 1
            if version_count == 0:
                FINAL_PATHOLOGIES[l] = 1
                DESIRED_PATHOLOGIES[l] = 1
            if version_count > 0:
                contents_file = f"{contents_dir}{file_prefix}{version_count}.contents"
                old_contents = dict()
                if os.path.isfile(contents_file):
                    print(f"Found an existing file for {alphabetical_l} run at: {contents_file}")
                    print(f"Checking for updates in {alphabetical_l} samples...")
                    with open(contents_file, 'r') as handle:
                        line_count = 0
                        for line in handle:
                            line_count += 1
                            line = line.strip("\n").split("\t")
                            if line_count == 1:
                                col_count = 0
                                for col in line:
                                    if col == "sample_id":
                                        id_col = col_count
                                    if col == "seq_type":
                                        seq_col = col_count
                                    col_count += 1
                            if line_count > 1:
                                sample_id = line[id_col]
                                seq_type = line[seq_col]
                                if not seq_type in old_contents:
                                    old_contents[seq_type] = [sample_id]
                                else:
                                    old_contents[seq_type].append(sample_id)

                new_contents = dict()
                for index, row in RUN_SAMPLES.iterrows():
                    sample_id = row["sample_id"]
                    seq_type = row["seq_type"]
                    pathology = row["pathology"]
                    if pathology in x and seq_type in SEQUENCIES:
                        if not seq_type in new_contents:
                            new_contents[seq_type] = [sample_id]
                        else:
                            new_contents[seq_type].append(sample_id)

                for sequence_type, new_samples in new_contents.items():
                    print("Looking for updates in sequence type...")
                    print(f"Checking if {sequence_type} in old contents...")
                    if not sequence_type in old_contents:
                        print(f"Current run list: {FINAL_PATHOLOGIES}")
                        FINAL_PATHOLOGIES[l] = version_count + 1
                        DESIRED_PATHOLOGIES[l] = version_count + 1
                        print(f"Found updates to seq_type: {sequence_type}!")
                        print(f"Run list is now: {FINAL_PATHOLOGIES}")
                        break # Break so it doesnt rerun again. If it reruns shouldn't affect the overalloutcome but at least it stops it from iterating more times?
                    elif sequence_type in old_contents:
                        print(f"Checking if any new sample updates in {sequence_type}...")
                        old_samples = old_contents[sequence_type]
                        if not all(sample in old_samples for sample in new_samples):
                            print("Found new sample ids!")
                            FINAL_PATHOLOGIES[l] = version_count + 1
                            DESIRED_PATHOLOGIES[l] = version_count + 1
                        elif all(sample in old_samples for sample in new_samples):
                            DESIRED_PATHOLOGIES[l] = version_count 

    assert len(FINAL_PATHOLOGIES) != 0, (
        f"Provided sample sets: {PATHOLOGIES} do not appear to be updated with new sample_ids or seq_types."
    )

    # Compare FINAL_PATHOLOGIES with DESIRED_PATHOLOGIES to see if still want to update
    # If FINAL_PATHOLOGY has NONE of the items in DESIRED_PATHOLOGY, means none are up to date and ignore

    not_updated = []
    for group in DESIRED_PATHOLOGIES:
        if group not in FINAL_PATHOLOGIES:
            not_updated.append(group)
    if len(not_updated) > 1:
        print(f"The following pathology/pathology groups have no updates and will not be run: {not_updated}")

    # Dictionary to pass all pathologies to R for merge script
    PATHOLOGY_DICT = {}
    for group in FINAL_PATHOLOGIES:
        PATHOLOGY_DICT[group] = {group: FINAL_PATHOLOGIES[group]}
    print(PATHOLOGY_DICT)

# Run with specific sample_ids
if CFG["options"]["sample_ids"] is not None:
    assert CFG["options"]["id_nickname"] is not None, (
        "If running specific sample_ids, please specify a nickname for your run!"
    )
    nickname = CFG["options"]["id_nickname"]
    CFG["options"]["mode"] = "samples"

    id_list = CFG["options"]["sample_ids"]

    # Should change PATHOLOGY_DICT to RUNS_DICT

    version_count = 0
    if os.path.exists(contents_dir):
        for f in os.listdir(contents_dir):
            contents_prefix = f"{nickname}_v"
            if contents_prefix in f and f.endswith(".contents"):
                version_count += 1
    
    if version_count > 0:
        print(f"Found {version_count} version(s) of run {nickname} at: {contents_prefix}{version_count}.contents\n Increasing the version number for this run, to make it easier to keep track ^-^")

    PATHOLOGY_DICT = {}
    PATHOLOGY_DICT[nickname] = {id_list: version_count + 1}

# Run with curated groups of samples subset
elif CFG["options"]["sample_subset"] is not None:
    subset_file = CFG["options"]["sample_subset_file"]
    subset = CFG["options"]["sample_subset_file"]
    CFG["options"]["mode"] = "subset"

    version_count = 0
    if os.path.exists(contents_dir):
        for f in os.listdir(contents_dir):
            contents_prefix = f"{subset}_v"
            if contents_prefix in f and f.endswith("contents"):
                version_count += 1
    
    PATHOLOGY_DICT = {}

    if version_count == 0:
        PATHOLOGY_DICT[subset] = {subset: 1}
    # Check contents to see if it has been updated?? 
    if version_count > 0:
        contents_file = f"{contents_dir}{contents_prefix}{version_count}.contents"
        old_contents = []
        with open(contents_file, 'r') as handle:
            line_counter = 0
            for line in handle:
                line = line.strip("\n").split("\t")
                line_counter += 1
                if line_counter == 1:
                    col_count = 0
                    for col in line:
                        if col == "sample_id":
                            sample_col = col_count
                        if col == subset:
                            subset_col = col_count
                        col_count += 1
                if line_counter > 1:
                    if line[subset_col] == 1:
                        old_contents.append(line[sample_col])
        new_contents = []
        with open(subset_file, 'r') as handle:
            line_counter = 0
            for line in handle:
                line = line.strip("\n").split("\t")
                line_counter += 1
                if line_counter == "1":
                    col_count = 0
                    for col in line:
                        if col == "sample_id":
                            sample_col = col_count
                        if col == subset:
                            subset_col = col_count
                        col_count += 1
                if line_counter > 1:
                    if line[subset_col] == "1":
                        new_contents.append(line[sample_col])
        
        if not all(samples in old_contents for samples in new_contents):
            PATHOLOGY_DICT[subset] = {subset: version_count + 1}

# Define rules to be run locally when using a compute cluster
localrules:
    _oncodriveclustl_get_maf,
    _run_oncodriveclustl,
    _oncodriveclustl_txt,
    _oncodriveclustl_all,


##### RULES #####

# Test that GAMBLR works:

md5hash = hashlib.md5()
conda_prefix = os.path.abspath(".snakemake/conda")
md5hash.update(conda_prefix.encode())
gamblr_file = open("/home/mcruz/lcr-modules/modules/oncodriveclustl/1.0.1/envs/GAMBLR.yaml", 'rb')
md5hash.update(gamblr_file.read())
gamblr_file.close()
h = md5hash.hexdigest()
gamblr_dir = conda_prefix + "/" + h[:8] + "/"

rule _install_gamblr:
    params:
        my_env = "/home/mcruz/miniconda3/envs/GAMBLR"
    output:
        installed = directory(gamblr_dir + "GAMBLR/lib/R/library/GAMBLR")
    shell:
        "cp -r {params.my_env} {output.installed}"

rule _oncodriveclustl_get_maf:
    input:
        # input fxn to look @ contents file, look for diffs, if no difference don't create content file / break run
        r_script = CFG["r_script"],
        gamblr = rules._install_gamblr.output.installed
    output:
        onco_input = CFG["dirs"]["inputs"] + "{genome_build}/{pathology}.CDS.maf"
    params:
        seq_type = CFG["options"]["seq_type"],
        genome_build = CFG["options"]["genome_build"], # can specify what build want to run
        sample_set = lambda wildcards: wildcards.pathology if config["lcr-modules"]["oncodriveclustl"]["options"]["mode"] == "pathology" else list(PATHOLOGY_DICT[wildcards.pathology]),
        options = CFG["options"]["clustl"],
        contents = lambda wildcards: expand(os.path.abspath(config["lcr-modules"]["oncodriveclustl"]["dirs"]["oncodriveclustl"] + wildcards.genome_build + "/contents/" + wildcards.pathology + "_v{version}.contents"), zip, version = PATHOLOGY_DICT[wildcards.pathology].values()) if config["lcr-modules"]["oncodriveclustl"]["options"]["mode"] in "pathology,samples,subset" else "",
        sample_metadata = "" if CFG["options"]["sample_metadata"] is None else "--metadata " + CFG["options"]["sample_metadata"],
        subset_file = "--setlist " + CFG["options"]["sample_subset_file"] if CFG["options"]["mode"] == "subset" else "",
        mode = CFG["options"]["mode"],
        output = os.path.abspath(CFG["dirs"]["inputs"] + "{genome_build}/{pathology}.CDS.maf")
    conda: CFG["conda_envs"]["r_env"]
    shell:
        op.as_one_line("""
            Rscript {input.r_script} 
            --sequencing "{params.seq_type}" --build {params.genome_build} --mode {params.mode}
            --sample_set "{params.sample_set}" --output "{params.output}" --contents {params.contents} {params.sample_metadata} {params.subset_file} {params.options}
        """)
    # script:
        # path to script 
-

rule _run_oncodriveclustl:
    input:
        # I might not need requirement for contents as input file since I have an assert above? Unless I want to manually check it in the shell/run command and proceed accordingly ://
        #contents = lambda wildcards: (config["lcr-modules"]["oncodriveclustl"]["dirs"]["inputs"] + "{genome_build}/{pathology}_v" + str(FINAL_PATHOLOGIES[wildcards.pathology]) + ".contents"),
        onco_input = rules._oncodriveclustl_get_maf.output.onco_input
    output:
        tsv = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{pathology}/clusters_results.tsv",
        txt = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{pathology}/elements_results.txt"
    log:
        stdout = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{pathology}/{pathology}_clustl_run.stdout.log",
        stderr = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{pathology}/{pathology}_clustl_run.stderr.log"
    conda: CFG["conda_envs"]["clustl"]
    params:
        local_path = "/projects/nhl_meta_analysis_scratch/gambl/ref_local/oncodrive/",
        output_dir = CFG["dirs"]["oncodriveclustl"] + "{genome_build}/{pathology}",
        threads = CFG["threads"]["clustl"],
        refs = CFG["options"]["reference"][CFG["options"]["genome_build"]],
        opts = CFG["options"]["clustl"]
    shell:
        op.as_one_line("""
       export BGDATA_LOCAL={params.local_path} && 
       oncodriveclustl -i {input.onco_input} -o {params.output_dir} -r {params.refs} --cores {params.threads} 
       {params.opts} > {log.stdout} 2> {log.stderr}
    """)

rule _oncodriveclustl_txt:
    input:
        tsv = rules._run_oncodriveclustl.output.tsv,
        txt = rules._run_oncodriveclustl.output.txt
    output:
        tsv = CFG["dirs"]["outputs"] + "{genome_build}/{pathology}.CDS.clusters.tsv",
        txt = CFG["dirs"]["outputs"] + "{genome_build}/{pathology}.CDS.elements.txt"
    run:
        op.relative_symlink(input.tsv, output.tsv)
        op.relative_symlink(input.txt, output.txt)

#rule _content_check:
 #   input:
#        CFG["dirs"]["outputs"] + "{{genome_build}}/{{pathology}}.CDS.clusters.tsv".format#(version = str(FINAL_PATHOLOGIES["{pathology}"]))
    #output:
    #    CFG["dirs"]["output"] + "{genome_build}/{pathology}.check"

#rule _oncodriveclustl_all:
#    input:
#        expand(
#            [
#                rules._oncodriveclustl_txt.output.tsv,
#                rules._oncodriveclustl_txt.output.txt
#            ],
#            zip,
#            pathology=list(PATHOLOGY_DICT),
#            genome_build=CFG["options"]["genome_build"]
#        )

rule _oncodriveclustl_all:
    input:
        expand(
            expand(
                [
                    rules._oncodriveclustl_txt.output.tsv,
                    rules._oncodriveclustl_txt.output.txt
                ],
                zip,
                genome_build=CFG["options"]["genome_build"],
                allow_missing=True
            ),
            pathology=list(PATHOLOGY_DICT)
        )
##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
