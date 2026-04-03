#!/usr/bin/env Rscript
#
# Notes: 
#   This script is intended for use with the dlbclone-1.0 module in LCR-modules. 
#   It uses the GAMBL.predict package associated metadata and binary feature matrix 
#   mutation data for training DLBCLone models.
#
#   This script runs UMAP (GAMBLR.predict::make_and_annotate_umap()) and optimizes 
#   parameters (GAMBLR.predict::DLBCLone_optimize_params) for classifying samples 
#   using a UMAP and k-nearest neighbors approach. A DLBCLone model object is then 
#   built and saved (GAMBLR.predict::DLBCLone_save_optimized) for future use in 
#   prediction. 
#
#   The user must appropriately set the config file with valid parameters if they 
#   choose to build a customized model, rather then the stored models already offered
#   by GAMBLR.predict. 

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressPackageStartupMessages({
  library(GAMBLR.predict)
  library(optparse)
  library(uwot)
})

# Determine arguments -----------------------------------------------------

# Define command line options
option_list <- list(
    make_option(c("-p", "--opt_model_path"), type = "character", help = "Output path for optimized model"),
    make_option(c("-n", "--model_name_prefix"), type = "character", help = "Name prefix for optimized model"),
    make_option(c("-b", "--core_features"), type = "character", help = "Core features as a list in R syntax"),
    make_option(c("-c", "--core_feature_multiplier"), type = "double", help = "Multiplier for core features during UMAP construction"),
    make_option(c("-d", "--hidden_features"), type = "character", help = "Hidden features as a vector in R syntax"),
    make_option(c("-e", "--truth_classes"), type = "character", help = "Truth classes as a vector in R syntax"),
    make_option(c("-g", "--min_k"), type = "integer", help = "Minimum k to test during optimization"),
    make_option(c("-i", "--max_k"), type = "integer", help = "Maximum k to test during optimization")
)

opt <- parse_args(OptionParser(option_list=option_list))

parse_nullable <- function(x) {
    if (x == "None") {
        return(NULL)
    }
    eval(parse(text = x))
}

opt_model_path <- opt$opt_model_path
model_name_prefix <- opt$model_name_prefix
core_features <- parse_nullable(opt$core_features)
core_feature_multiplier <- opt$core_feature_multiplier
hidden_features <- parse_nullable(opt$hidden_features)
truth_classes <- eval(parse(text=opt$truth_classes))
min_k <- opt$min_k
max_k <- opt$max_k

# Load GAMBL.predict metadata and binary feature matrix--------------------
metadata <- readr::read_tsv(system.file("extdata/dlbcl_meta_with_dlbclass.tsv",package = "GAMBLR.predict")) %>%
    dplyr::filter(lymphgen %in% truth_classes)

mutation_matrix <- readr::read_tsv(system.file("extdata/all_full_status.tsv",package = "GAMBLR.predict")) %>%
  tibble::column_to_rownames("sample_id") 

# Run GAMBLR.predict::make_and_annotate_umap()-----------------------------
umap_model <- make_and_annotate_umap(
    df = mutation_matrix,
    metadata = metadata,
    core_features = core_features,
    core_feature_multiplier = core_feature_multiplier,
    hidden_features = hidden_features
)

# Run GAMBLR.predict::DLBCLone_optimize_params()----------------------------
opt_model <- DLBCLone_optimize_params(  
    combined_mutation_status_df = umap_model$features, 
    metadata_df = metadata,
    umap_out = umap_model,
    truth_classes = truth_classes,
    min_k = min_k,
    max_k = max_k
)

# Output --------------------------------------------------------------------------------------

# saving optimized model
DLBCLone_save_optimized( 
    DLBCLone_model = opt_model,
    base_path = opt_model_path,
    name_prefix = model_name_prefix
)
