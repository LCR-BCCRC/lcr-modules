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
    subdirectories = ["inputs", "dlbclone", "outputs"],
)

# TODO: Parallelization. What if you are predicting on multiple models? What if you are building multiple models?
# TODO: Appropriately incorporate str() as required by LCR modules
# TODO: Appropriately assign CFG "dirs" + "outputs"
# TODO: Unit tests for assemble_genetic_features rule and its incorporation into dlbclone_predict rule


##### RULES #####
if "launch_date" in CFG:
    launch_date = CFG['launch_date']
else:
    launch_date = datetime.today().strftime('%Y-%m')


#  Convert CFG Python dictionary lists to an R list of c(...) vectors
def as_r_list(d):
    if d is None:
        return "NULL"
    items = []
    for k, v in d.items():
        r_vec = "c(" + ",".join(f'"{x}"' for x in v) + ")"
        items.append(f"{k}={r_vec}")
    return "list(" + ",".join(items) + ")"

# Convert CFG Python dictionary to an R named list
def as_r_named_list(d):
    if d is None:
        return "NULL"
    items = []
    for k, v in d.items():
        items.append(f'{k}="{v}"')
    return "list(" + ",".join(items) + ")"

# Convert CFG Python lists to an R c(...) vector
def as_r_c(d):
    if d is None or len(d) == 0:
        return "None"
    return "c(" + ",".join(f'"{x}"' for x in d) + ")"


# Define rules to be run locally when using a compute cluster
localrules:
    _dlbclone_sif_pull,
    #_dlbclone_input_maf,
    _dlbclone_build_model,
    _dlbclone_assemble_genetic_features,
    _dlbclone_predict,
    _dlbclone_all

# Pull SIF file from lcr-sifs
rule _dlbclone_sif_pull:
    output:
        sif = CFG["dirs"]["dlbclone"] + "gamblr.sif" #dlbclone.sif
    #log:
     #   stdout = CFG["logs"]["dlbclone"] + "install.stdout.log",
      #  stderr = CFG["logs"]["dlbclone"] + "install.stderr.log"
 #   conda:
  #      CFG["conda"]["singularity"]
    shell:
        """
        singularity pull {output.sif} docker://lklossok/gamblr:latest
        """


# STEP 1: INPUT SYMLINKS
# Symlinks the input files into the module results directory (under '00-inputs/')
#rule _dlbclone_input_maf:
 #   input:
  #      maf = CFG["inputs"]["sample_maf"] 
   # output:
    #    maf = CFG["dirs"]["inputs"] + "maf/{seq_type}--{genome_build}/{tumour_id}--{normal_id}--{pair_status}.maf"
#    group:
 #       "lymphgen"
  #  run:
   #     op.relative_symlink(input.maf, output.maf)


# Uses only GAMBL metadata and binary mutation data as training for DLBCLone models
rule _dlbclone_build_model: 
    input:
        dlbclone_model = CFG["inputs"]["dlbclone_model"], 
        sif = str(rules._dlbclone_sif_pull.output.sif)
    output:
        opt_model_rds_path = CFG["inputs"]["opt_model_path"] + "/" + CFG["inputs"]["model_name_prefix"] + "_model.rds",
        opt_model_uwot_path = CFG["inputs"]["opt_model_path"] + "/" + CFG["inputs"]["model_name_prefix"] + "_umap.uwot"
    params:
        opt_model_path = CFG["inputs"]["opt_model_path"],
        model_name_prefix = CFG["inputs"]["model_name_prefix"],
        core_features = as_r_list(CFG["options"]["core_features"]),
        core_feature_multiplier = CFG["options"]["core_feature_multiplier"],
        hidden_features = as_r_c(CFG["options"]["hidden_features"]),
        truth_classes = as_r_c(CFG["options"]["truth_classes"]),
        min_k = CFG["options"]["min_k"],
        max_k = CFG["options"]["max_k"]
    #log:
     #   CFG["logs"]["opt_model_path"] + "/" + CFG["inputs"]["model_name_prefix"] + "_model.rds",
      #  CFG["logs"]["opt_model_path"] + "/" + CFG["inputs"]["model_name_prefix"] + "_umap.uwot"
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    shell:
        op.as_one_line("""
            apptainer exec --cleanenv {input.sif} Rscript {input.dlbclone_model}  
                --opt_model_path {params.opt_model_path} 
                --model_name_prefix {params.model_name_prefix}
                --core_features '{params.core_features}' 
                --core_feature_multiplier {params.core_feature_multiplier} 
                --hidden_features '{params.hidden_features}' 
                --truth_classes '{params.truth_classes}' 
                --min_k {params.min_k} 
                --max_k {params.max_k}
        """)


# Determine if custom dlbclone model(s) need to be generated for predictions
USE_GAMBLR = CFG["inputs"]["opt_model_path"] == "GAMBLR.predict"

HAS_CUSTOM_MODELS = False
if not USE_GAMBLR:
    model_dir = Path(CFG["inputs"]["opt_model_path"])
    HAS_CUSTOM_MODELS = any(model_dir.glob("*_model.rds"))

if USE_GAMBLR:
    print("Using GAMBLR.predict pre-made dlbclone models")
elif HAS_CUSTOM_MODELS:
    print("Using custom pre-made dlbclone model(s)")
else:
    print("Building custom dlbclone model(s)")

def model_inputs():
    if USE_GAMBLR or HAS_CUSTOM_MODELS:
        return []
    else:
        return rules._dlbclone_build_model.output


rule _dlbclone_assemble_genetic_features:
    input:
        test_metadata = CFG["inputs"]["test_metadata_dir"],
        test_maf = CFG["inputs"]["test_maf_dir"],
        dlbclone_assemble_genetic_features = CFG["inputs"]["dlbclone_assemble_genetic_features"], 
        sif = str(rules._dlbclone_sif_pull.output.sif)
    output:
        test_mutation_matrix = CFG["inputs"]["test_matrix_dir"] + "/" + CFG["inputs"]["matrix_name"] + "_mutation_matrix.tsv"
    params:
        sv_from_metadata = as_r_named_list(CFG["options"]["sv_from_metadata"]),
        translocation_status = as_r_named_list(CFG["options"]["translocation_status"]),
        maf_sample_id_colname = CFG["options"]["maf_sample_id_colname"],
        metadata_sample_id_colname = CFG["options"]["metadata_sample_id_colname"],
        truth_column_colname = CFG["options"]["truth_column_colname"],
        genes = as_r_c(CFG["options"]["genes"]),
        synon_genes = as_r_c(CFG["options"]["synon_genes"]),
        hotspot_genes = as_r_c(CFG["options"]["hotspot_genes"]),
        sv_value = CFG["options"]["sv_value"],
        synon_value = CFG["options"]["synon_value"],
        coding_value = CFG["options"]["coding_value"],
        include_ashm = CFG["options"]["include_ashm"],
        annotated_sv = CFG["options"]["annotated_sv"],
        include_GAMBL_sv = CFG["options"]["include_GAMBL_sv"]
    #log:
    #   CFG["logs"]["dlbclone_predict"]["test_data_dir"]
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    shell:
        op.as_one_line("""
        apptainer exec --cleanenv {input.sif} Rscript {input.dlbclone_assemble_genetic_features} 
            --metadata {input.test_metadata} 
            --maf {input.test_maf} 
            --maf_sample_id_colname {params.maf_sample_id_colname}
            --metadata_sample_id_colname {params.metadata_sample_id_colname} 
            --sv_from_metadata '{params.sv_from_metadata}' 
            --translocation_status '{params.translocation_status}' 
            --truth_column_colname {params.truth_column_colname}
            --output_matrix_dir {output.test_mutation_matrix} 
            --genes '{params.genes}'
            --synon_genes '{params.synon_genes}'
            --hotspot_genes '{params.hotspot_genes}'
            --sv_value {params.sv_value}
            --synon_value {params.synon_value}
            --coding_value {params.coding_value}
            --include_ashm {params.include_ashm}
            --annotated_sv {params.annotated_sv}
            --include_GAMBL_sv {params.include_GAMBL_sv}
        """)


rule _dlbclone_predict:
    input:
        models = model_inputs(), 
        mutation_matrix = str(rules._dlbclone_assemble_genetic_features.output.test_mutation_matrix), #CFG["inputs"]["test_matrix_dir"], # first column is sample ID and all other columns are features
        dlbclone_predict = CFG["inputs"]["dlbclone_predict"], 
        sif = str(rules._dlbclone_sif_pull.output.sif)
    output:
        predictions = CFG["inputs"]["pred_dir"] + "/" + CFG["inputs"]["model_name_prefix"] + "--{launch_date}/{matrix_name}_DLBCLone_predictions.tsv"
    params:
        opt_model_path = CFG["inputs"]["opt_model_path"],
        model_name_prefix = CFG["inputs"]["model_name_prefix"],
        fill_missing = CFG["options"]["fill_missing"], 
        drop_extra = CFG["options"]["drop_extra"]
    #log:
     #   CFG["logs"]["pred_dir"] + "/" + CFG["inputs"]["model_name_prefix"] + "--{launch_date}/{matrix_name}_DLBCLone_predictions.tsv"
    #container:
     #   "docker://lklossok/dlbclone:latest"
    threads:
        CFG["threads"]["quant"]
    resources:
        **CFG["resources"]["quant"]
    shell:
        op.as_one_line("""
        apptainer exec --cleanenv {input.sif} Rscript {input.dlbclone_predict} 
            --test_features {input.mutation_matrix} 
            --output_dir {output.predictions} 
            --model_path {params.opt_model_path} 
            --model_prefix {params.model_name_prefix} 
            --fill_missing {params.fill_missing} 
            --drop_extra {params.drop_extra}
        """)


rule _dlbclone_all:
    input:
        expand(
            [
                rules._dlbclone_assemble_genetic_features.output.test_mutation_matrix #_dlbclone_predict.output.predictions
            ],
            #model_name_prefix = CFG["inputs"]["model_name_prefix"],
            matrix_name = CFG["inputs"]["matrix_name"],
            launch_date = launch_date
        )

##### CLEANUP #####


# Perform some clean-up tasks, including storing the module-specific
# configuration on disk and deleting the `CFG` variable
op.cleanup_module(CFG)
