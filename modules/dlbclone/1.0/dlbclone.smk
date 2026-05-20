#!/usr/bin/env snakemake


##### ATTRIBUTION #####


# Original Author:  Luke Klossok
# Module Author:    Luke Klossok
# Contributors:     N/A


##### SETUP #####

# Import package with useful functions for developing analysis modules
import oncopipe as op
from datetime import datetime
from pathlib import Path
from glob import glob

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
# `CFG` is a shortcut to `config["lcr-modules"]["dlbclone"]`
CFG = op.setup_module(
    name = "dlbclone",
    version = "1.0",
    # TODO: If applicable, add more granular output subdirectories
    subdirectories = ["inputs", "prepare_inputs", "dlbclone", "outputs"],
)


# Determine if custom dlbclone model(s) need to be generated for predictions
USE_GAMBLR = CFG["inputs"]["opt_model_path"] == "GAMBLR.predict"

if not USE_GAMBLR:
    MODEL_DIR = Path(CFG["inputs"]["opt_model_path"])
    MODEL_PREFIX = CFG["inputs"]["model_name_prefix"]
    MODEL_RDS = str(MODEL_DIR / f"{MODEL_PREFIX}_model.rds")
    MODEL_UWOT = str(MODEL_DIR / f"{MODEL_PREFIX}_umap.uwot")


# Define rules to be run locally when using a compute cluster
if USE_GAMBLR:
    localrules:
        _dlbclone_input_metadata,
        _dlbclone_input_maf,
        _dlbclone_input_sv,
        get_combined_sv,
        curate_sv,
        combine_maf,
        combine_sv_metadata,
        _dlbclone_assemble_genetic_features,
        _dlbclone_build_model,
        _dlbclone_predict,
        _dlbclone_all


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')

#  Convert CFG Python dictionary lists to an R list of c(...) vectors
def as_r_list(d):
    if d is None:
        return "None"
    items = []
    for k, v in d.items():
        r_vec = "c(" + ",".join(f'"{x}"' for x in v) + ")"
        items.append(f"{k}={r_vec}")
    return "list(" + ",".join(items) + ")"

# Convert CFG Python dictionary to an R named list
def as_r_named_list(d):
    if d is None:
        return "None"
    items = []
    for k, v in d.items():
        items.append(f'{k}="{v}"')
    return "list(" + ",".join(items) + ")"

# Convert CFG Python lists to an R c(...) vector
def as_r_c(d):
    if d is None or len(d) == 0:
        return "None"
    return "c(" + ",".join(f'"{x}"' for x in d) + ")"


# STEP 1: INPUT SYMLINKS
# Symlinks the input files into the module results directory (under '00-inputs/')

rule _dlbclone_input_metadata:
    input:
        metadata = CFG["inputs"]["test_metadata_dir"] 
    output:
        metadata = CFG["dirs"]["inputs"] + "metadata/{sample_set}--{launch_date}_metadata.tsv"
    run:
        op.relative_symlink(input.metadata, output.metadata)

rule _dlbclone_input_maf:
    input:
        maf = CFG["inputs"]["sample_maf"] 
    output:
        maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
    run:
        op.relative_symlink(input.maf, output.maf)

rule _dlbclone_input_sv: 
    input: 
        sv = CFG["inputs"]["sample_sv"]  if CFG["inputs"]["sample_sv"] != "" else CFG["inputs"]["empty_sv"]
    output: 
        sv = CFG["dirs"]["inputs"] + "sv/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.bedpe"
    run: 
        op.relative_symlink(input.sv, output.sv)

# STEP 2: PREPARE AGGREGATE INPUTS

rule combine_maf:
    input:
        mafs = expand(
            CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf",
            zip,
            seq_type=SAMPLES.seq_type,
            genome_build=SAMPLES.genome_build,
            tumour_id=SAMPLES.tumour_id,
            normal_id=SAMPLES.normal_id,
            pair_status=SAMPLES.pair_status
        )
    output:
        maf = temp(CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/combined.maf")
    shell:
        op.as_one_line("""
            head -n 1 {input.mafs[0]} > {output.maf} &&
            tail -n +2 -q {input.mafs} >> {output.maf}
        """)

rule get_combined_sv:
    input:
        svs = expand(
            rules._dlbclone_input_sv.output.sv,
            zip,
            seq_type=SAMPLES.seq_type,
            genome_build=SAMPLES.genome_build,
            tumour_id=SAMPLES.tumour_id,
            normal_id=SAMPLES.normal_id,
            pair_status=SAMPLES.pair_status
        ) if CFG["inputs"]["sample_combined_sv"] == "" else [],
        combined = CFG["inputs"]["sample_combined_sv"] if CFG["inputs"]["sample_combined_sv"] != "" else []
    output:
        sv = temp(CFG["dirs"]["inputs"] + "sv/{sample_set}--{launch_date}/combined_input_sv.bedpe")
    run:
        if input.combined:
            op.relative_symlink(input.combined, output.sv)
        else:
            shell("""
                head -n 1 {input.svs[0]} > {output.sv} &&
                tail -n +2 -q {input.svs} >> {output.sv}
            """)


# STEP 3 PREPARE SV & FEATURE MATRIX INPUTS

rule curate_sv:
    input:
        script=CFG["scripts"]["prepare_svs"], 
        sv=str(rules.get_combined_sv.output.sv),
        metadata=str(rules._dlbclone_input_metadata.output.metadata)
    output:
        output_metadata=CFG["dirs"]["prepare_inputs"] + "{sample_set}--{launch_date}/combined_metadata.tsv"
    params:
        real_bcl2=as_r_c(CFG["options"]["real_bcl2"]),
        real_bcl6=as_r_c(CFG["options"]["real_bcl6"]),
        metadata_FISH_cols=as_r_named_list(CFG["inputs"]["metadata_FISH_cols"])
    log:
        CFG["logs"]["outputs"] + "{sample_set}--{launch_date}/combined_metadata.tsv"
    singularity: 
        "docker://lklossok/predict:latest"
    shell:
        op.as_one_line("""
            Rscript {input.script}
                --sv_input {input.sv}
                --metadata {input.metadata}
                --output_metadata {output.output_metadata}
                --real_bcl2 '{params.real_bcl2}' 
                --real_bcl6 '{params.real_bcl6}'
                --sv_fish_colname '{params.metadata_FISH_cols}'
            2> {log}
        """)

rule _dlbclone_assemble_genetic_features:
    input:
        script = CFG["scripts"]["dlbclone_assemble_genetic_features"],
        test_metadata =  str(rules.curate_sv.output.output_metadata), 
        test_maf = str(rules.combine_maf.output.maf) 
    output:
        test_mutation_matrix = CFG["dirs"]["outputs"] + "{sample_set}--{launch_date}/{matrix_name}_mutation_matrix.tsv"
    params:
        sv_from_metadata = as_r_named_list(CFG["options"]["sv_from_metadata"]),
        translocation_status = as_r_named_list(CFG["options"]["translocation_status"]),
        maf_sample_id_colname = CFG["options"]["maf_sample_id_colname"],
        metadata_sample_id_colname = CFG["options"]["metadata_sample_id_colname"],
        genes_coding = as_r_c(CFG["options"]["genes_coding"]),
        genes_noncoding = as_r_c(CFG["options"]["genes_noncoding"]),
        genes_hotspot = as_r_c(CFG["options"]["genes_hotspot"]),
        genes_driver = as_r_c(CFG["options"]["genes_driver"]),
        sv_value = CFG["options"]["sv_value"],
        noncoding_value = CFG["options"]["noncoding_value"],
        coding_value = CFG["options"]["coding_value"],
        driver_value = CFG["options"]["driver_value"]
    log:
        CFG["logs"]["outputs"] + "{sample_set}--{launch_date}/{matrix_name}_mutation_matrix.log"
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    singularity: 
            "docker://lklossok/predict:latest"
    shell:
        op.as_one_line("""
            Rscript {input.script} 
                --metadata {input.test_metadata} 
                --maf {input.test_maf} 
                --maf_sample_id_colname {params.maf_sample_id_colname}
                --metadata_sample_id_colname {params.metadata_sample_id_colname} 
                --sv_from_metadata '{params.sv_from_metadata}' 
                --translocation_status '{params.translocation_status}'
                --output_matrix_dir {output.test_mutation_matrix} 
                --genes_coding '{params.genes_coding}'
                --genes_noncoding '{params.genes_noncoding}'
                --genes_hotspot '{params.genes_hotspot}'
                --genes_driver '{params.genes_driver}'
                --sv_value {params.sv_value}
                --noncoding_value {params.noncoding_value}
                --coding_value {params.coding_value}
                --driver_value {params.driver_value}
            2> {log}
        """)


# STEP 4: RUN DLBCLONE RULES 

if not USE_GAMBLR:
    # Uses only GAMBL metadata and binary mutation data as training for DLBCLone models
    rule _dlbclone_build_model: 
        input:
            script = CFG["scripts"]["dlbclone_model"]
        output:
            model_rds = MODEL_RDS,
            model_uwot = MODEL_UWOT
        params:
            opt_model_path = CFG["inputs"]["opt_model_path"],
            model_name_prefix = CFG["inputs"]["model_name_prefix"],
            core_features = as_r_list(CFG["options"]["core_features"]),
            core_feature_multiplier = CFG["options"]["core_feature_multiplier"],
            hidden_features = as_r_c(CFG["options"]["hidden_features"]),
            truth_classes = as_r_c(CFG["options"]["truth_classes"]),
            min_k = CFG["options"]["min_k"],
            max_k = CFG["options"]["max_k"]
        threads:
            CFG["threads"]["quant"]
        resources:
            **CFG["resources"]["quant"]
        singularity: 
            "docker://lklossok/predict:latest"
        shell:
            op.as_one_line("""
                Rscript {input.script}  
                    --opt_model_path {params.opt_model_path} 
                    --model_name_prefix {params.model_name_prefix}
                    --core_features '{params.core_features}' 
                    --core_feature_multiplier {params.core_feature_multiplier} 
                    --hidden_features '{params.hidden_features}' 
                    --truth_classes '{params.truth_classes}' 
                    --min_k {params.min_k} 
                    --max_k {params.max_k}
            """)


rule _dlbclone_predict:
    input:
        script = CFG["scripts"]["dlbclone_predict"],
        models = [] if USE_GAMBLR else [MODEL_RDS, MODEL_UWOT],
        mutation_matrix = str(rules._dlbclone_assemble_genetic_features.output.test_mutation_matrix) 
    output:
        predictions = CFG["dirs"]["outputs"] + "{sample_set}--{launch_date}/{model_name_prefix}--{matrix_name}_DLBCLone_predictions.tsv"
    params:
        opt_model_path = CFG["inputs"]["opt_model_path"],
        model_name_prefix = CFG["inputs"]["model_name_prefix"],
        fill_missing = CFG["options"]["fill_missing"], 
        drop_extra = CFG["options"]["drop_extra"],
    log:
        CFG["logs"]["outputs"] + "{sample_set}--{launch_date}/{model_name_prefix}--{matrix_name}_DLBCLone_predictions.tsv"
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    singularity: 
            "docker://lklossok/predict:latest"
    shell:
        op.as_one_line("""
            Rscript {input.script} 
                --test_features {input.mutation_matrix} 
                --output_dir {output.predictions}
                --model_path {params.opt_model_path} 
                --model_prefix {params.model_name_prefix} 
                --fill_missing {params.fill_missing} 
                --drop_extra {params.drop_extra}
            2> {log}
        """)


rule _dlbclone_all:
    input:
        expand(
            [
                str(rules.curate_sv.output.output_metadata) #str(rules._dlbclone_predict.output.predictions)
            ],
            zip,
            sample_id = SAMPLES["sample_id"],
            patient_id = SAMPLES["patient_id"],
            tissue_status = SAMPLES["tissue_status"],
            seq_type=SAMPLES["seq_type"],
            genome_build=SAMPLES["genome_build"],
            tumour_id=SAMPLES["tumour_id"],
            normal_id=SAMPLES["normal_id"],
            pair_status=SAMPLES["pair_status"], 
            allow_missing = True,
            launch_date = launch_date,
            sample_set = CFG["inputs"]["sample_set"],
            matrix_name = CFG["inputs"]["matrix_name"],
            model_name_prefix = CFG["inputs"]["model_name_prefix"]
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
